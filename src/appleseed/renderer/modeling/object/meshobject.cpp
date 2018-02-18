
//
// This source file is part of appleseed.
// Visit http://appleseedhq.net/ for additional information and resources.
//
// This software is released under the MIT license.
//
// Copyright (c) 2010-2013 Francois Beaune, Jupiter Jazz Limited
// Copyright (c) 2014-2017 Francois Beaune, The appleseedhq Organization
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
#include "meshobject.h"

// appleseed.renderer headers.
#include "renderer/kernel/tessellation/statictessellation.h"
#include "renderer/modeling/object/iregion.h"
#include "renderer/modeling/object/meshobjectprimitives.h"
#include "renderer/modeling/object/meshobjectreader.h"
#include "renderer/modeling/object/triangle.h"

// appleseed.foundation headers.
#include "foundation/utility/api/apiarray.h"
#include "foundation/utility/api/apistring.h"
#include "foundation/utility/api/specializedapiarrays.h"
#include "foundation/utility/containers/dictionary.h"
#include "foundation/utility/foreach.h"
#include "foundation/utility/string.h"

// Standard headers.
#include <cassert>
#include <string>
#include <vector>

using namespace foundation;
using namespace std;

namespace renderer
{

//
// MeshObject class implementation.
//

namespace
{
    const char* Model = "mesh_object";

    // A region that simply wraps a static tessellation.
    class TessWrappingRegion
      : public IRegion
    {
      public:
        explicit TessWrappingRegion(StaticTriangleTess* tess)
          : m_tess(tess)
          , m_lazy_tess(tess)
        {
        }

        GAABB3 compute_local_bbox() const override
        {
            return m_tess->compute_local_bbox();
        }

        Lazy<StaticTriangleTess>& get_static_triangle_tess() const override
        {
            return m_lazy_tess;
        }

      private:
        StaticTriangleTess*                 m_tess;
        mutable Lazy<StaticTriangleTess>    m_lazy_tess;
    };
}

struct MeshObject::Impl
{
    bool                        m_is_displaced;

    StaticTriangleTess          m_tess;
    TessWrappingRegion          m_region;
    RegionKit                   m_region_kit;
    mutable Lazy<RegionKit>     m_lazy_region_kit;

    StaticTriangleTess          m_displaced_tess;
    TessWrappingRegion          m_displaced_region;
    RegionKit                   m_displaced_region_kit;
    mutable Lazy<RegionKit>     m_lazy_displaced_region_kit;

    vector<string>              m_material_slots;

    Impl()
      : m_is_displaced(false)
      , m_region(&m_tess)
      , m_lazy_region_kit(&m_region_kit)
      , m_displaced_region(&m_displaced_tess)
      , m_lazy_displaced_region_kit(&m_displaced_region_kit)
    {
        m_region_kit.push_back(&m_region);
        m_displaced_region_kit.push_back(&m_displaced_region);
    }

    static GScalar compute_height(const GVector2& uv)
    {
        //return (1.0f + sin(norm(GVector2(0.5f) - uv) * 2.0f)) * 0.005f;
        return uv[0] * 0.05f;
    }

    void build_displaced_surface(const Source* map, const Source* multiplier)
    {
        // Build a direct relationship between vertices, vertex normals and texture coordinates.
        vector<GVector3> vertex_normals(m_tess.m_vertices.size(), GVector3(0.0f));
        vector<GVector2> vertex_uvs(m_tess.m_vertices.size(), GVector2(0.0f));
        for (const Triangle& tri : m_tess.m_primitives)
        {
            if (tri.m_n0 != Triangle::None) vertex_normals[tri.m_v0] = m_tess.m_vertex_normals[tri.m_n0];
            if (tri.m_n1 != Triangle::None) vertex_normals[tri.m_v1] = m_tess.m_vertex_normals[tri.m_n1];
            if (tri.m_n2 != Triangle::None) vertex_normals[tri.m_v2] = m_tess.m_vertex_normals[tri.m_n2];
            if (tri.m_a0 != Triangle::None) vertex_uvs[tri.m_v0] = m_tess.get_tex_coords(tri.m_a0);
            if (tri.m_a1 != Triangle::None) vertex_uvs[tri.m_v1] = m_tess.get_tex_coords(tri.m_a1);
            if (tri.m_a2 != Triangle::None) vertex_uvs[tri.m_v2] = m_tess.get_tex_coords(tri.m_a2);
        }

        // Displace and push the original vertices into the final tessellation.
        for (size_t i = 0, e = m_tess.m_vertices.size(); i < e; ++i)
        {
            const GVector3& n = vertex_normals[i];
            const GVector2& uv = vertex_uvs[i];
            GVector3 v = m_tess.m_vertices[i];

            // Apply displacement.
            assert(is_normalized(n));
            v += n * compute_height(uv);

            m_displaced_tess.m_vertices.push_back(v);
        }

        // Push the original normals into the final tessellation.
        for (const GVector3& n : m_tess.m_vertex_normals)
            m_displaced_tess.m_vertex_normals.push_back(n);

        // Split and push the original triangles into the final tessellation.
        for (const Triangle& tri : m_tess.m_primitives)
            split_or_insert_triangle(tri, vertex_normals, vertex_uvs, 0);
    }

    void split_or_insert_triangle(
        const Triangle&     tri,
        vector<GVector3>&   vertex_normals,
        vector<GVector2>&   vertex_uvs,
        const size_t        level)
    {
        const GVector3& v0 = m_displaced_tess.m_vertices[tri.m_v0];
        const GVector3& v1 = m_displaced_tess.m_vertices[tri.m_v1];
        const GVector3& v2 = m_displaced_tess.m_vertices[tri.m_v2];

        const GScalar square_max_edge_length = square(0.02f);

        const GScalar sn01 = square_norm(v1 - v0);
        const GScalar sn02 = square_norm(v2 - v0);
        const GScalar sn12 = square_norm(v2 - v1);

        const bool split =
            sn01 > square_max_edge_length ||
            sn02 > square_max_edge_length ||
            sn12 > square_max_edge_length;

        if (split)
        {
            // Vertex indices of the edge centers.
            const size_t i01 = m_displaced_tess.m_vertices.size() + 0;
            const size_t i02 = m_displaced_tess.m_vertices.size() + 1;
            const size_t i12 = m_displaced_tess.m_vertices.size() + 2;

            // Compute and push normals of edge centers.
            const GVector3& n0 = vertex_normals[tri.m_v0];
            const GVector3& n1 = vertex_normals[tri.m_v1];
            const GVector3& n2 = vertex_normals[tri.m_v2];
            const GVector3 n01 = normalize(n0 + n1);  // todo: improve?
            const GVector3 n02 = normalize(n0 + n2);  // todo: improve?
            const GVector3 n12 = normalize(n1 + n2);  // todo: improve?
            ensure_minimum_size(vertex_normals, i12 + 1);
            vertex_normals[i01] = n01;
            vertex_normals[i02] = n02;
            vertex_normals[i12] = n12;

            // Compute and push UV coordinates of edge centers.
            const GVector2& uv0 = vertex_uvs[tri.m_v0];
            const GVector2& uv1 = vertex_uvs[tri.m_v1];
            const GVector2& uv2 = vertex_uvs[tri.m_v2];
            const GVector2 uv01 = 0.5f * (uv0 + uv1);
            const GVector2 uv02 = 0.5f * (uv0 + uv2);
            const GVector2 uv12 = 0.5f * (uv1 + uv2);
            ensure_minimum_size(vertex_uvs, i12 + 1);
            vertex_uvs[i01] = uv01;
            vertex_uvs[i02] = uv02;
            vertex_uvs[i12] = uv12;

            // Compute initial position of edge centers.
            GVector3 v01 = 0.5f * (v0 + v1);
            GVector3 v02 = 0.5f * (v0 + v2);
            GVector3 v12 = 0.5f * (v1 + v2);

            // Apply displacement to edge centers.
            v01 += n01 * compute_height(uv01);
            v02 += n02 * compute_height(uv02);
            v12 += n12 * compute_height(uv12);

            // Push edge centers.
            m_displaced_tess.m_vertices.push_back(v01);
            m_displaced_tess.m_vertices.push_back(v02);
            m_displaced_tess.m_vertices.push_back(v12);

            if (level < 8)
            {
                // Recursively consider the four new triangles.
                split_or_insert_triangle(Triangle(tri.m_v0, i01, i02), vertex_normals, vertex_uvs, level + 1);
                split_or_insert_triangle(Triangle(tri.m_v1, i12, i01), vertex_normals, vertex_uvs, level + 1);
                split_or_insert_triangle(Triangle(tri.m_v2, i02, i12), vertex_normals, vertex_uvs, level + 1);
                split_or_insert_triangle(Triangle(i01, i12, i02),      vertex_normals, vertex_uvs, level + 1);
            }
        }
        else
        {
            m_displaced_tess.m_primitives.push_back(tri);
        }
    }
};

MeshObject::MeshObject(
    const char*             name,
    const ParamArray&       params)
  : Object(name, params)
  , impl(new Impl())
{
    m_inputs.declare("alpha_map", InputFormatFloat, "");
    m_inputs.declare("displacement_map", InputFormatFloat, "");
    m_inputs.declare("displacement_multiplier", InputFormatFloat, "1.0");
}

MeshObject::~MeshObject()
{
    delete impl;
}

void MeshObject::release()
{
    delete this;
}

const char* MeshObject::get_model() const
{
    return Model;
}

bool MeshObject::on_render_begin(
    const Project&          project,
    IAbortSwitch*           abort_switch)
{
    if (!Object::on_render_begin(project, abort_switch))
        return false;

    impl->m_is_displaced = false;

    const Source* displacement_map = m_inputs.source("displacement_map");
    if (displacement_map != nullptr)
    {
        if (get_tex_coords_count() > 0)
        {
            impl->build_displaced_surface(
                displacement_map,
                m_inputs.source("displacement_multiplier"));
            impl->m_is_displaced = true;
        }
        else
        {
            RENDERER_LOG_WARNING(
                "mesh object \"%s\" does not have texture coordinates, cannot apply displacement.",
                get_path().c_str());
        }
    }

    return true;
}

void MeshObject::on_render_end(const Project& project)
{
    if (impl->m_is_displaced)
    {
        impl->m_displaced_tess.clear_release_memory();
        impl->m_is_displaced = false;
    }

    Object::on_render_end(project);
}

bool MeshObject::on_frame_begin(
    const Project&          project,
    const BaseGroup*        parent,
    OnFrameBeginRecorder&   recorder,
    IAbortSwitch*           abort_switch)
{
    if (!Object::on_frame_begin(project, parent, recorder, abort_switch))
        return false;

    m_alpha_map = get_uncached_alpha_map();

    return true;
}

void MeshObject::on_frame_end(
    const Project&          project,
    const BaseGroup*        parent)
{
    m_alpha_map = nullptr;
    Object::on_frame_end(project, parent);
}

bool MeshObject::has_alpha_map() const
{
    if (!m_params.strings().exist("alpha_map"))
        return false;

    const char* value = m_params.strings().get("alpha_map");

    return !is_empty_string(value);
}

const Source* MeshObject::get_uncached_alpha_map() const
{
    return m_inputs.source("alpha_map");
}

GAABB3 MeshObject::compute_local_bbox() const
{
    return impl->m_tess.compute_local_bbox();
}

Lazy<RegionKit>& MeshObject::get_region_kit()
{
    return impl->m_is_displaced
        ? impl->m_lazy_displaced_region_kit
        : impl->m_lazy_region_kit;
}

void MeshObject::reserve_vertices(const size_t count)
{
    impl->m_tess.m_vertices.reserve(count);
}

size_t MeshObject::push_vertex(const GVector3& vertex)
{
    const size_t index = impl->m_tess.m_vertices.size();
    impl->m_tess.m_vertices.push_back(vertex);
    return index;
}

size_t MeshObject::get_vertex_count() const
{
    return impl->m_tess.m_vertices.size();
}

const GVector3& MeshObject::get_vertex(const size_t index) const
{
    return impl->m_tess.m_vertices[index];
}

void MeshObject::reserve_vertex_normals(const size_t count)
{
    impl->m_tess.m_vertex_normals.reserve(count);
}

size_t MeshObject::push_vertex_normal(const GVector3& normal)
{
    assert(is_normalized(normal));

    const size_t index = impl->m_tess.m_vertex_normals.size();
    impl->m_tess.m_vertex_normals.push_back(normal);
    return index;
}

size_t MeshObject::get_vertex_normal_count() const
{
    return impl->m_tess.m_vertex_normals.size();
}

const GVector3& MeshObject::get_vertex_normal(const size_t index) const
{
    return impl->m_tess.m_vertex_normals[index];
}

void MeshObject::clear_vertex_normals()
{
    impl->m_tess.m_vertex_normals.clear();
}

void MeshObject::reserve_vertex_tangents(const size_t count)
{
    impl->m_tess.reserve_vertex_tangents(count);
}

size_t MeshObject::push_vertex_tangent(const GVector3& tangent)
{
    return impl->m_tess.push_vertex_tangent(tangent);
}

size_t MeshObject::get_vertex_tangent_count() const
{
    return impl->m_tess.get_vertex_tangent_count();
}

GVector3 MeshObject::get_vertex_tangent(const size_t index) const
{
    return impl->m_tess.get_vertex_tangent(index);
}

void MeshObject::reserve_tex_coords(const size_t count)
{
    impl->m_tess.reserve_tex_coords(count);
}

size_t MeshObject::push_tex_coords(const GVector2& tex_coords)
{
    return impl->m_tess.push_tex_coords(tex_coords);
}

size_t MeshObject::get_tex_coords_count() const
{
    return impl->m_tess.get_tex_coords_count();
}

GVector2 MeshObject::get_tex_coords(const size_t index) const
{
    return impl->m_tess.get_tex_coords(index);
}

void MeshObject::reserve_triangles(const size_t count)
{
    impl->m_tess.m_primitives.reserve(count);
}

size_t MeshObject::push_triangle(const Triangle& triangle)
{
    const size_t index = impl->m_tess.m_primitives.size();
    impl->m_tess.m_primitives.push_back(triangle);
    return index;
}

size_t MeshObject::get_triangle_count() const
{
    return impl->m_tess.m_primitives.size();
}

const Triangle& MeshObject::get_triangle(const size_t index) const
{
    return impl->m_tess.m_primitives[index];
}

Triangle& MeshObject::get_triangle(const size_t index)
{
    return impl->m_tess.m_primitives[index];
}

void MeshObject::clear_triangles()
{
    impl->m_tess.m_primitives.clear();
}

void MeshObject::set_motion_segment_count(const size_t count)
{
    impl->m_tess.set_motion_segment_count(count);
}

size_t MeshObject::get_motion_segment_count() const
{
    return impl->m_tess.get_motion_segment_count();
}

void MeshObject::set_vertex_pose(
    const size_t            vertex_index,
    const size_t            motion_segment_index,
    const GVector3&         vertex)
{
    impl->m_tess.set_vertex_pose(vertex_index, motion_segment_index, vertex);
}

GVector3 MeshObject::get_vertex_pose(
    const size_t            vertex_index,
    const size_t            motion_segment_index) const
{
    return impl->m_tess.get_vertex_pose(vertex_index, motion_segment_index);
}

void MeshObject::clear_vertex_poses()
{
    impl->m_tess.clear_vertex_poses();
}

void MeshObject::set_vertex_normal_pose(
    const size_t            normal_index,
    const size_t            motion_segment_index,
    const GVector3&         normal)
{
    impl->m_tess.set_vertex_normal_pose(normal_index, motion_segment_index, normal);
}

GVector3 MeshObject::get_vertex_normal_pose(
    const size_t            normal_index,
    const size_t            motion_segment_index) const
{
    return impl->m_tess.get_vertex_normal_pose(normal_index, motion_segment_index);
}

void MeshObject::clear_vertex_normal_poses()
{
    impl->m_tess.clear_vertex_normal_poses();
}

void MeshObject::set_vertex_tangent_pose(
    const size_t            tangent_index,
    const size_t            motion_segment_index,
    const GVector3&         tangent)
{
    impl->m_tess.set_vertex_tangent_pose(tangent_index, motion_segment_index, tangent);
}

GVector3 MeshObject::get_vertex_tangent_pose(
    const size_t            tangent_index,
    const size_t            motion_segment_index) const
{
    return impl->m_tess.get_vertex_tangent_pose(tangent_index, motion_segment_index);
}

void MeshObject::clear_vertex_tangent_poses()
{
    impl->m_tess.clear_vertex_tangent_poses();
}

void MeshObject::reserve_material_slots(const size_t count)
{
    impl->m_material_slots.reserve(count);
}

size_t MeshObject::push_material_slot(const char* name)
{
    const size_t index = impl->m_material_slots.size();
    impl->m_material_slots.emplace_back(name);
    return index;
}

size_t MeshObject::get_material_slot_count() const
{
    return impl->m_material_slots.size();
}

const char* MeshObject::get_material_slot(const size_t index) const
{
    return impl->m_material_slots[index].c_str();
}

void MeshObject::collect_asset_paths(StringArray& paths) const
{
    if (m_params.strings().exist("filename"))
        paths.push_back(m_params.get("filename"));
    else if (m_params.dictionaries().exist("filename"))
    {
        const StringDictionary& filepaths = m_params.dictionaries().get("filename").strings();
        for (const_each<StringDictionary> i = filepaths; i; ++i)
            paths.push_back(i->value());
    }
}

void MeshObject::update_asset_paths(const StringDictionary& mappings)
{
    if (m_params.strings().exist("filename"))
        m_params.set("filename", mappings.get(m_params.get("filename")));
    else if (m_params.dictionaries().exist("filename"))
    {
        StringDictionary& filepaths = m_params.dictionaries().get("filename").strings();
        for (const_each<StringDictionary> i = filepaths; i; ++i)
            filepaths.set(i->key(), mappings.get(i->value()));
    }
}


//
// MeshObjectFactory class implementation.
//

void MeshObjectFactory::release()
{
    delete this;
}

const char* MeshObjectFactory::get_model() const
{
    return Model;
}

Dictionary MeshObjectFactory::get_model_metadata() const
{
    return
        Dictionary()
            .insert("name", Model)
            .insert("label", "Mesh Object");
}

DictionaryArray MeshObjectFactory::get_input_metadata() const
{
    DictionaryArray metadata;

    metadata.push_back(
        Dictionary()
            .insert("name", "alpha_map")
            .insert("label", "Alpha Map")
            .insert("type", "colormap")
            .insert("entity_types",
                Dictionary()
                    .insert("color", "Colors")
                    .insert("texture_instance", "Textures"))
            .insert("min",
                Dictionary()
                    .insert("value", "0.0")
                    .insert("type", "hard"))
            .insert("max",
                Dictionary()
                    .insert("value", "1.0")
                    .insert("type", "hard"))
            .insert("use", "optional"));

    metadata.push_back(
        Dictionary()
            .insert("name", "displacement_map")
            .insert("label", "Displacement Map")
            .insert("type", "colormap")
            .insert("entity_types",
                Dictionary()
                    .insert("texture_instance", "Textures"))
            .insert("min",
                Dictionary()
                    .insert("value", "-1.0")
                    .insert("type", "soft"))
            .insert("max",
                Dictionary()
                    .insert("value", "1.0")
                    .insert("type", "soft"))
            .insert("use", "optional"));

    metadata.push_back(
        Dictionary()
            .insert("name", "displacement_multiplier")
            .insert("label", "Displacement Multiplier")
            .insert("type", "numeric")
            .insert("min",
                Dictionary()
                    .insert("value", "0.0")
                    .insert("type", "soft"))
            .insert("max",
                Dictionary()
                    .insert("value", "1.0")
                    .insert("type", "soft"))
            .insert("use", "optional")
            .insert("default", "1.0"));

    return metadata;
}

auto_release_ptr<Object> MeshObjectFactory::create(
    const char*             name,
    const ParamArray&       params) const
{
    return auto_release_ptr<Object>(new MeshObject(name, params));
}

bool MeshObjectFactory::create(
    const char*             name,
    const ParamArray&       params,
    const SearchPaths&      search_paths,
    const bool              omit_loading_assets,
    ObjectArray&            objects) const
{
    if (params.strings().exist("primitive"))
    {
        auto_release_ptr<MeshObject> mesh = create_primitive_mesh(name, params);
        if (mesh.get() == nullptr)
            return false;

        objects.push_back(mesh.release());
        return true;
    }

    if (omit_loading_assets)
    {
        objects.push_back(create(name, params).release());
        return true;
    }

    MeshObjectArray object_array;
    if (!MeshObjectReader::read(
            search_paths,
            name,
            params,
            object_array))
        return false;

    objects = array_vector<ObjectArray>(object_array);
    return true;
}

}   // namespace renderer
