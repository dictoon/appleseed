
//
// This source file is part of appleseed.
// Visit https://appleseedhq.net/ for additional information and resources.
//
// This software is released under the MIT license.
//
// Copyright (c) 2010-2013 Francois Beaune, Jupiter Jazz Limited
// Copyright (c) 2014-2018 Francois Beaune, The appleseedhq Organization
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//

// Interface header.
#include "backwardlightsampler.h"

// appleseed.renderer headers.
#include "renderer/global/globallogger.h"
#include "renderer/kernel/intersection/intersector.h"
#include "renderer/kernel/lighting/directlightingintegrator.h"
#include "renderer/kernel/lighting/lightsample.h"
#include "renderer/kernel/lighting/pathtracer.h"
#include "renderer/kernel/lighting/pathvertex.h"
#include "renderer/kernel/lighting/scatteringmode.h"
#include "renderer/kernel/lighting/tracer.h"
#include "renderer/kernel/shading/shadingpoint.h"
#include "renderer/kernel/texturing/texturecache.h"
#include "renderer/modeling/camera/camera.h"
#include "renderer/modeling/edf/edf.h"
#include "renderer/modeling/frame/frame.h"
#include "renderer/modeling/light/light.h"
#include "renderer/modeling/material/material.h"
#include "renderer/modeling/scene/scene.h"

// appleseed.foundation headers.
#include "foundation/image/canvasproperties.h"
#include "foundation/image/image.h"
#include "foundation/math/knn.h"
#include "foundation/math/mis.h"
#include "foundation/math/rng/distribution.h"
#include "foundation/platform/compiler.h"
#include "foundation/platform/defaulttimers.h"
#include "foundation/utility/makevector.h"
#include "foundation/utility/memory.h"
#include "foundation/utility/otherwise.h"
#include "foundation/utility/poison.h"

// Standard headers.
#include <cassert>
#include <string>

using namespace foundation;
using namespace std;

namespace renderer
{

//
// BackwardLightSampler class implementation.
//

Dictionary BackwardLightSampler::get_params_metadata()
{
    Dictionary metadata;

    metadata.insert(
        "algorithm",
        Dictionary()
            .insert("type", "enum")
            .insert("values", "cdf|lighttree|lightprobes")
            .insert("default", "cdf")
            .insert("label", "Light Sampler")
            .insert("help", "Light sampling algoritm")
            .insert(
                "options",
                Dictionary()
                    .insert(
                        "cdf",
                        Dictionary()
                            .insert("label", "CDF")
                            .insert("help", "Cumulative Distribution Function"))
                    .insert(
                        "lighttree",
                        Dictionary()
                            .insert("label", "Light Tree")
                            .insert("help", "Lights organized in a BVH"))
                    .insert(
                        "lightprobes",
                        Dictionary()
                            .insert("label", "Light Probes")
                            .insert("help", "Precomputed Light Probes"))));

    metadata.merge(LightSamplerBase::get_params_metadata());

    return metadata;
}

BackwardLightSampler::BackwardLightSampler(
    const Scene&                        scene,
    const ParamArray&                   params)
  : LightSamplerBase(params)
  , m_scene(scene)
{
    // Retrieve the sampling algorithm.
    const string algorithm_str =
        params.get_optional<string>("algorithm", "cdf", make_vector("cdf", "lighttree", "lightprobes"));
    m_algorithm =
        algorithm_str == "cdf" ? Algorithm::CDF :
        algorithm_str == "lighttree" ? Algorithm::LightTree :
        Algorithm::LightProbes;

    RENDERER_LOG_INFO("collecting light emitters...");

    // Collect all non-physical lights and separate them according to their
    // compatibility with the LightTree.
    collect_non_physical_lights(
        scene.assembly_instances(),
        TransformSequence(),
        [&](const NonPhysicalLightInfo& light_info)
        {
            if (m_algorithm == Algorithm::LightTree
                && ((light_info.m_light->get_flags() & Light::LightTreeCompatible) != 0))
            {
                // Insert into light tree compatible lights.
                m_light_tree_lights.push_back(light_info);
            }
            else
            {
                // Insert into non-physical lights to be evaluated using a CDF.
                const size_t light_index = m_non_physical_lights.size();
                m_non_physical_lights.push_back(light_info);

                // Insert the light into the CDF.
                // todo: compute importance.
                float importance = 1.0f;
                importance *= light_info.m_light->get_uncached_importance_multiplier();
                m_non_physical_lights_cdf.insert(light_index, importance);
            }
        });
    m_non_physical_light_count = m_non_physical_lights.size();

    // Collect all light-emitting shapes.
    collect_emitting_shapes(
        scene.assembly_instances(),
        TransformSequence(),
        [&](
            const Material* material,
            const float     area,
            const size_t    emitting_shape_index)
        {
            if (m_algorithm == Algorithm::LightTree)
            {
                // Only accept this shape if its material has an EDF.
                // This excludes shapes with light-emitting OSL materials
                // since these are not handled by the light tree yet.
                return material->get_uncached_edf() != nullptr;
            }
            else
            {
                // Retrieve the EDF and get the importance multiplier.
                float importance_multiplier = 1.0f;
                if (const EDF* edf = material->get_uncached_edf())
                    importance_multiplier = edf->get_uncached_importance_multiplier();

                // Compute the probability density of this shape.
                const float shape_importance = m_params.m_importance_sampling ? static_cast<float>(area) : 1.0f;
                const float shape_prob = shape_importance * importance_multiplier;

                // Insert the light-emitting shape into the CDF.
                m_emitting_shapes_cdf.insert(emitting_shape_index, shape_prob);

                // Accept this shape.
                return true;
            }
        });

    // Build the hash table of emitting shapes.
    build_emitting_shape_hash_table();

    // Prepare the non-physical lights CDF for sampling.
    if (m_non_physical_lights_cdf.valid())
        m_non_physical_lights_cdf.prepare();

    if (m_algorithm == Algorithm::LightTree)
    {
        // Initialize the LightTree only after the lights are collected.
        m_light_tree.reset(new LightTree(m_light_tree_lights, m_emitting_shapes));

        // Build the light tree.
        const vector<size_t> tri_index_to_node_index = m_light_tree->build();
        assert(tri_index_to_node_index.size() == m_emitting_shapes.size());

        // Associate light tree nodes to emitting shapes.
        for (size_t i = 0, e = m_emitting_shapes.size(); i < e; ++i)
            m_emitting_shapes[i].m_light_tree_node_index = tri_index_to_node_index[i];
    }
    else
    {
        // Prepare the light-emitting shapes CDF for smapling.
        if (m_emitting_shapes_cdf.valid())
            m_emitting_shapes_cdf.prepare();

        // Store the shape probability densities into the emitting shapes.
        for (size_t i = 0, e = m_emitting_shapes.size(); i < e; ++i)
            m_emitting_shapes[i].m_shape_prob = m_emitting_shapes_cdf[i].second;
    }

    RENDERER_LOG_INFO(
        "found %s %s, %s %s, %s emitting %s.",
        pretty_int(m_non_physical_light_count).c_str(),
        plural(m_non_physical_light_count, "non-physical light").c_str(),
        pretty_int(m_light_tree_lights.size() + m_emitting_shapes.size()).c_str(),
        plural(m_light_tree_lights.size() + m_emitting_shapes.size(), "light-tree compatible light").c_str(),
        pretty_int(m_emitting_shapes.size()).c_str(),
        plural(m_emitting_shapes.size(), "shape").c_str());
}

namespace
{
    // Direct lighting sampling probes settings.
    constexpr size_t MaxShapeCount = 64;
    constexpr size_t PixelStride = 16;
    constexpr size_t LightSampleCount = 1;
    constexpr size_t NeighborProbeCount = 1;
    constexpr float ResidualWeightFraction = 0.05f;
}

class BackwardLightSampler::PathVisitor
{
  public:
    using CDFType = CDF<size_t, float>;
    using ContribType = std::pair<size_t, float>;

    BackwardLightSampler*           m_light_sampler;
    const ShadingContext&           m_shading_context;
    SamplingContext&                m_sampling_context;
    const bool                      m_enable_caustics;
    const size_t                    m_light_sample_count;
    vector<Vector3d>&               m_light_probe_positions;
    vector<CDFType>&                m_emitting_shape_cdfs;
    bool                            m_is_indirect_lighting;
    vector<ContribType>             m_unoccluded_contribs;
    vector<float>                   m_weights;

    PathVisitor(
        BackwardLightSampler*       light_sampler,
        const ShadingContext&       shading_context,
        SamplingContext&            sampling_context,
        const bool                  enable_caustics,
        const size_t                light_sample_count,
        vector<Vector3d>&           light_probe_positions,
        vector<CDFType>&            emitting_shape_cdfs)
      : m_light_sampler(light_sampler)
      , m_shading_context(shading_context)
      , m_sampling_context(sampling_context)
      , m_enable_caustics(enable_caustics)
      , m_light_sample_count(light_sample_count)
      , m_light_probe_positions(light_probe_positions)
      , m_emitting_shape_cdfs(emitting_shape_cdfs)
      , m_is_indirect_lighting(false)
    {
    }

    void on_first_diffuse_bounce(const PathVertex& vertex)
    {
    }
    
    bool accept_scattering(
        const ScatteringMode::Mode  prev_mode,
        const ScatteringMode::Mode  next_mode) const
    {
        assert(next_mode != ScatteringMode::None);

        if (!m_enable_caustics)
        {
            // Don't follow paths leading to caustics.
            if (ScatteringMode::has_diffuse_or_volume(prev_mode) &&
                ScatteringMode::has_glossy_or_specular(next_mode))
                return false;
        }

        return true;
    }
    
    void on_miss(const PathVertex& vertex)
    {
    }

    void on_hit(const PathVertex& vertex)
    {
        assert(vertex.m_shading_point != nullptr);

        if (vertex.m_bsdf == nullptr)
            return;

        CDFType cdf;

        const BSDFSampler bsdf_sampler(
            *vertex.m_bsdf,
            vertex.m_bsdf_data,
            0,                          // BSDF sampling modes -- not used by the methods we use
            *vertex.m_shading_point);

        const DirectLightingIntegrator dl_integrator(
            m_shading_context,
            *m_light_sampler,
            bsdf_sampler,
            vertex.get_time(),
            0,                          // light sampling modes -- not used by the methods we use
            0,                          // material sample count -- not used by the methods we use
            m_light_sample_count,       // not used by the methods we use
            0.0f,                       // not used by the methods we use
            m_is_indirect_lighting);

        const size_t shape_count = m_light_sampler->m_emitting_shapes.size();
        cdf.reserve(min(shape_count, MaxShapeCount));

        if (shape_count <= MaxShapeCount)
        {
            for (size_t shape_index = 0; shape_index < shape_count; ++shape_index)
            {
                m_shading_context.get_arena().clear();
                m_sampling_context.split_in_place(2, m_light_sample_count);

                Spectrum irradiance;
                irradiance.set(0.0f);

                for (size_t i = 0, e = m_light_sample_count; i < e; ++i)
                {
                    LightSample light_sample;
                    light_sample.m_light = nullptr;
                    m_light_sampler->sample_emitting_shape(
                        vertex.get_time(),
                        m_sampling_context.next2<Vector2f>(),
                        shape_index,
                        1.0f,               // shape probability
                        light_sample);

                    dl_integrator.add_emitting_shape_sample_contribution_to_irradiance(
                        light_sample,
                        true,               // occluded
                        irradiance);
                }

                cdf.insert(shape_index, sum_value(irradiance));
            }
        }
        else
        {
            // Compute unoccluded contribution for all shapes.
            m_unoccluded_contribs.reserve(shape_count);
            for (size_t shape_index = 0; shape_index < shape_count; ++shape_index)
            {
                m_shading_context.get_arena().clear();
                m_sampling_context.split_in_place(2, m_light_sample_count);

                Spectrum irradiance;
                irradiance.set(0.0f);

                for (size_t i = 0, e = m_light_sample_count; i < e; ++i)
                {
                    LightSample light_sample;
                    light_sample.m_light = nullptr;
                    m_light_sampler->sample_emitting_shape(
                        vertex.get_time(),
                        m_sampling_context.next2<Vector2f>(),
                        shape_index,
                        1.0f,               // shape probability
                        light_sample);

                    dl_integrator.add_emitting_shape_sample_contribution_to_irradiance(
                        light_sample,
                        false,              // unoccluded
                        irradiance);
                }

                m_unoccluded_contribs.emplace_back(shape_index, sum_value(irradiance));
            }

            // Sort shapes by decreasing contribution.
            sort(
                m_unoccluded_contribs.begin(),
                m_unoccluded_contribs.end(),
                [](const ContribType& lhs, const ContribType& rhs)
                {
                    return rhs.second < lhs.second;
                });

            float top_weight = 0.0f;
            for (size_t i = 0; i < MaxShapeCount; ++i)
                top_weight += m_unoccluded_contribs[i].second;

            const size_t discarded_shape_count = shape_count - MaxShapeCount;
            const float residual_weight = top_weight * ResidualWeightFraction;
            const float discarded_shape_weight = residual_weight / discarded_shape_count;
            for (size_t i = MaxShapeCount; i < shape_count; ++i)
                m_unoccluded_contribs[i].second = discarded_shape_weight;

            for (const ContribType& contrib : m_unoccluded_contribs)
                cdf.insert(contrib.first, contrib.second);

            /*ensure_minimum_size(m_weights, shape_count);

            for (size_t shape_index = 0; shape_index < shape_count; ++shape_index)
                m_weights[shape_index] = 0.0f;

            for (size_t shape_sample_index = 0; shape_sample_index < MaxShapeCount; ++shape_sample_index)
            {
                const size_t shape_index = static_cast<size_t>(rand_int1(m_rng, 0, static_cast<int32>(shape_count - 1)));

                m_shading_context.get_arena().clear();

                m_sampling_context.split_in_place(2, m_light_sample_count);

                Spectrum irradiance;
                irradiance.set(0.0f);

                for (size_t i = 0, e = m_light_sample_count; i < e; ++i)
                {
                    LightSample light_sample;
                    light_sample.m_light = nullptr;
                    m_light_sampler->sample_emitting_shape(
                        vertex.get_time(),
                        m_sampling_context.next2<Vector2f>(),
                        shape_index,
                        1.0f,               // shape probability
                        light_sample);

                    dl_integrator.add_emitting_shape_sample_contribution_to_irradiance(
                        light_sample,
                        irradiance);
                }

                m_weights[shape_index] = sum_value(irradiance);
            }

            float total_weight = 0.0f;
            size_t discarded_shape_count = 0;

            for (size_t shape_index = 0; shape_index < shape_count; ++shape_index)
            {
                const float shape_weight = m_weights[shape_index];
                total_weight += shape_weight;
                if (shape_weight == 0.0f)
                    ++discarded_shape_count;
            }

            assert(shape_count > MaxShapeCount);
            assert(discarded_shape_count >= shape_count - MaxShapeCount);

            const float residual_weight = total_weight > 0.0f ? ResidualWeightFraction * total_weight : 1.0f;
            const float discarded_shape_weight = residual_weight / discarded_shape_count;

            for (size_t shape_index = 0; shape_index < shape_count; ++shape_index)
            {
                const float shape_weight = m_weights[shape_index];
                cdf.insert(shape_index, shape_weight > 0.0f ? shape_weight : discarded_shape_weight);
            }*/
        }

        assert(!cdf.empty());

        if (cdf.valid())
        {
            cdf.prepare();
            m_emitting_shape_cdfs.push_back(cdf);

            m_light_probe_positions.push_back(vertex.m_shading_point->get_point());
        }
    }

    void on_scatter(PathVertex& vertex)
    {
        assert(vertex.m_scattering_modes != ScatteringMode::None);

        // Any light contribution after a diffuse or glossy bounce is considered indirect.
        if (ScatteringMode::has_diffuse_or_glossy_or_volume(vertex.m_prev_mode))
            m_is_indirect_lighting = true;

        // When caustics are disabled, disable glossy and specular components after a diffuse or volume bounce.
        // Note that accept_scattering() is later going to return false in this case.
        if (!m_enable_caustics)
        {
            if (vertex.m_prev_mode == ScatteringMode::Diffuse ||
                vertex.m_prev_mode == ScatteringMode::Volume)
                vertex.m_scattering_modes &= ~(ScatteringMode::Glossy | ScatteringMode::Specular);
        }
    }
};

namespace
{
    struct VolumeVisitor
    {
        bool accept_scattering(
            const ScatteringMode::Mode  prev_mode)
        {
            return true;
        }
    
        void on_scatter(PathVertex& vertex)
        {
        }

        void visit_ray(PathVertex& vertex, const ShadingRay& volume_ray)
        {
        }
    };
}

void BackwardLightSampler::create_sampling_probes(
    const Frame&                        frame,
    const TraceContext&                 trace_context,
    TextureStore&                       texture_store,
    OIIOTextureSystem&                  oiio_texture_system,
    OSLShadingSystem&                   osl_shading_system,
    IAbortSwitch&                       abort_switch)
{
    assert(m_emitting_shape_cdfs.empty());

    if (m_algorithm != Algorithm::LightProbes)
        return;

    // No light source in the scene.
    if (!has_lights())  
        return;

    vector<Vector3d> light_probe_positions;

    TextureCache texture_cache(texture_store);
    const Intersector intersector(trace_context, texture_cache);
    Arena arena;
    OSLShaderGroupExec shadergroup_exec(osl_shading_system, arena);
    Tracer tracer(m_scene, intersector, shadergroup_exec);

    const ShadingContext shading_context(
        intersector,
        tracer,
        texture_cache,
        oiio_texture_system,
        shadergroup_exec,
        arena,
        0);     // thread index

    const Image& image = frame.image();
    const CanvasProperties& props = image.properties();

    const size_t location_count =
        truncate<size_t>(ceil(static_cast<double>(props.m_canvas_width) / PixelStride)) *
        truncate<size_t>(ceil(static_cast<double>(props.m_canvas_height) / PixelStride));

    RENDERER_LOG_INFO(
        "creating direct lighting sampling probes at %s location%s for %s light-emitting shape%s...",
        pretty_uint(location_count).c_str(),
        location_count > 1 ? "s" : "",
        pretty_uint(m_emitting_shapes.size()).c_str(),
        m_emitting_shapes.size() > 1 ? "s" : "");

    size_t location_done = 0;
    size_t last_reported_completed_percent = 0;

    for (size_t y = 0; y < props.m_canvas_height; y += PixelStride)
    {
        for (size_t x = 0; x < props.m_canvas_width; x += PixelStride)
        {
            if (abort_switch.is_aborted())
                return;

            // Create a sampling context.
            const size_t pixel_index = y * props.m_canvas_width + x;
            const size_t instance = hash_uint32(static_cast<uint32>(pixel_index));
            SamplingContext::RNGType rng(instance, instance);
            SamplingContext sampling_context(
                rng,
                SamplingContext::QMCMode,
                instance);                  // initial instance number

            // Compute the sample position in NDC.
            const Vector2d sample_position = frame.get_sample_position(x + 0.5, y + 0.5);

            // 1/4 of a pixel, like in RenderMan RIS.
            const Vector2d sample_position_dx(1.0 / (4.0 * props.m_canvas_width), 0.0);
            const Vector2d sample_position_dy(0.0, -1.0 / (4.0 * props.m_canvas_height));

            // Construct a primary ray.
            ShadingRay primary_ray;
            m_scene.get_active_camera()->spawn_ray(
                sampling_context,
                Dual2d(sample_position, sample_position_dx, sample_position_dy),
                primary_ray);

            // Trace the ray.
            ShadingPoint shading_point;
            intersector.trace(primary_ray, shading_point);

            // todo: get from the path tracer config?
            const bool EnableCaustics = false;
            const size_t RRMinPathLength = 1000;
            const size_t MaxBounces = 3;
            const size_t MaxDiffuseBounces = 3;
            const size_t MaxGlossyBounces = 3;
            const size_t MaxSpecularBounces = 3;
            const size_t MaxVolumeBounces = 3;
            const bool ClampRoughness = false;

            PathVisitor path_visitor(
                this,
                shading_context,
                sampling_context,
                EnableCaustics,
                LightSampleCount,
                light_probe_positions,
                m_emitting_shape_cdfs);

            VolumeVisitor volume_visitor;

            PathTracer<PathVisitor, VolumeVisitor, false> path_tracer(
                path_visitor,
                volume_visitor,
                RRMinPathLength,
                MaxBounces,
                MaxDiffuseBounces,
                MaxGlossyBounces,
                MaxSpecularBounces,
                MaxVolumeBounces,
                ClampRoughness);

            path_tracer.trace(
                sampling_context,
                shading_context,
                shading_point);

            ++location_done;

            const size_t completed_percent = location_done * 100 / location_count;
            if (completed_percent - last_reported_completed_percent >= 5)
            {
                last_reported_completed_percent = completed_percent;
                RENDERER_LOG_INFO(
                    "creating direct lighting sampling probes, %s completed...",
                    pretty_percent(location_done, location_count).c_str());
            }
        }
    }

    assert(location_done == location_count);
    assert(light_probe_positions.size() == m_emitting_shape_cdfs.size());

    RENDERER_LOG_INFO(
        "created %s direct lighting sampling probe%s.",
        pretty_uint(light_probe_positions.size()).c_str(),
        light_probe_positions.size() > 1 ? "s" : "");

    if (!light_probe_positions.empty())
    {
        // todo: use Builder::build_move_points()?
        knn::Builder3d builder(m_light_probe_tree);
        builder.build<DefaultWallclockTimer>(
            &light_probe_positions[0],
            light_probe_positions.size());
    }
}

void BackwardLightSampler::sample_lightset(
    const ShadingRay::Time&             time,
    const Vector3f&                     s,
    const ShadingPoint&                 shading_point,
    LightSample&                        light_sample) const
{
    switch (m_algorithm)
    {
      case Algorithm::CDF:
        sample_emitting_shapes(time, s, light_sample);
        break;

      case Algorithm::LightTree:
        sample_light_tree(time, s, shading_point, light_sample);
        break;

      case Algorithm::LightProbes:
        sample_light_probes(time, s, shading_point, light_sample);
        break;

      assert_otherwise;
    }
}

float BackwardLightSampler::evaluate_pdf(
    const ShadingPoint&                 light_shading_point,
    const ShadingPoint&                 surface_shading_point) const
{
    assert(light_shading_point.is_triangle_primitive());

    const EmittingShapeKey shape_key(
        light_shading_point.get_assembly_instance().get_uid(),
        light_shading_point.get_object_instance_index(),
        light_shading_point.get_primitive_index());

    const auto* shape_ptr = m_emitting_shape_hash_table.get(shape_key);

    if (shape_ptr == nullptr)
        return 0.0f;

    const EmittingShape* shape = *shape_ptr;

    float shape_prob;
    poison(shape_prob);

    switch (m_algorithm)
    {
      case Algorithm::CDF:
        shape_prob = shape->m_shape_prob;
        break;

      case Algorithm::LightTree:
        shape_prob =
            m_light_tree->evaluate_node_pdf(
                surface_shading_point,
                shape->m_light_tree_node_index);
        break;

      case Algorithm::LightProbes:
        shape_prob = evaluate_light_probes_pdf(shape, surface_shading_point);
        break;

      assert_otherwise;
    }

    const float pdf = shape_prob * shape->m_rcp_area;
    assert(pdf >= 0.0f);

    return pdf;
}

void BackwardLightSampler::sample_light_tree(
    const ShadingRay::Time&             time,
    const Vector3f&                     s,
    const ShadingPoint&                 shading_point,
    LightSample&                        light_sample) const
{
    assert(has_lightset());

    LightType light_type;
    size_t light_index;
    float light_prob;
    m_light_tree->sample(
        shading_point,
        s[0],
        light_type,
        light_index,
        light_prob);

    if (light_type == NonPhysicalLightType)
    {
        // Fetch the light.
        const NonPhysicalLightInfo& light_info = m_light_tree_lights[light_index];
        light_sample.m_light = light_info.m_light;

        // Evaluate and store the transform of the light.
        light_sample.m_light_transform =
              light_info.m_light->get_transform()
            * light_info.m_transform_sequence.evaluate(time.m_absolute);

        // Store the probability density of this light.
        light_sample.m_probability = light_prob;
    }
    else
    {
        assert(light_type == EmittingShapeType);
        sample_emitting_shape(
            time,
            Vector2f(s[1], s[2]),
            light_index,
            light_prob,
            light_sample);
    }

    assert(light_sample.m_light || light_sample.m_shape);
    assert(light_sample.m_probability > 0.0f);
}

void BackwardLightSampler::sample_light_probes(
    const ShadingRay::Time&             time,
    const Vector3f&                     s,
    const ShadingPoint&                 shading_point,
    LightSample&                        light_sample) const
{
    if (s[0] < 0.5f)
    {
        sample_emitting_shapes(
            time,
            Vector3f(2.0f * s[0], s[1], s[2]),
            light_sample);
    }
    else
    {
        // todo: this allocates memory!
        knn::Answer<double> answer(NeighborProbeCount);

        knn::Query3d query(m_light_probe_tree, answer);
        query.run(shading_point.get_point());
        assert(answer.size() == 1);

        // answer.sort();

        const size_t probe_index = m_light_probe_tree.remap(answer.get(0).m_index);
        assert(probe_index < m_emitting_shape_cdfs.size());

        const CDFType& emitting_shapes_cdf = m_emitting_shape_cdfs[probe_index];
        assert(emitting_shapes_cdf.valid());

        const EmitterCDF::ItemWeightPair result = emitting_shapes_cdf.sample(2.0f * (s[0] - 0.5f));
        const size_t shape_index = result.first;
        const float shape_prob = result.second;

        light_sample.m_light = nullptr;
        sample_emitting_shape(
            time,
            Vector2f(s[1], s[2]),
            shape_index,
            shape_prob,
            light_sample);

        assert(light_sample.m_shape);
        assert(light_sample.m_probability > 0.0f);
    }
}

float BackwardLightSampler::evaluate_light_probes_pdf(
    const EmittingShape*                shape,
    const ShadingPoint&                 surface_shading_point) const
{
    // todo: this allocates memory!
    knn::Answer<double> answer(NeighborProbeCount);

    knn::Query3d query(m_light_probe_tree, answer);
    query.run(surface_shading_point.get_point());
    assert(answer.size() == 1);

    // answer.sort();

    const size_t probe_index = m_light_probe_tree.remap(answer.get(0).m_index);
    assert(probe_index < m_emitting_shape_cdfs.size());

    const CDFType& emitting_shapes_cdf = m_emitting_shape_cdfs[probe_index];
    assert(emitting_shapes_cdf.valid());

    for (size_t i = 0, e = emitting_shapes_cdf.size(); i < e; ++i)
    {
        const CDFType::ItemWeightPair& item = emitting_shapes_cdf[i];
        if (&m_emitting_shapes[item.first] == shape)
        {
            // MIS, single sample model, balance heuristic.
            return 0.5f * item.second + 0.5f * shape->m_shape_prob;
        }
    }

    APPLESEED_UNREACHABLE;
    return 0.0f;
}

}   // namespace renderer
