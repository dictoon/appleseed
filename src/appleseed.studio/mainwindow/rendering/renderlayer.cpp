
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
#include "renderlayer.h"

// appleseed.renderer headers.
#include "renderer/global/globallogger.h"

// appleseed.foundation headers.
#include "foundation/image/canvasproperties.h"
#include "foundation/image/image.h"
#include "utility/gl.h"

// Qt headers.
#include <QMutexLocker>
#include <QOpenGLFunctions_4_1_Core>
#include <QOpenGLTexture>

using namespace foundation;

namespace appleseed {
namespace studio {

//
// RenderLayer class implementation.
//

RenderLayer::RenderLayer(
    const std::size_t       width,
    const std::size_t       height,
    OCIO::ConstConfigRcPtr  ocio_config,
    QWidget*                parent)
  : RenderWidget(width, height, ocio_config, parent)
  , m_gl_initialized(false)
  , m_refresh_gl_texture(false)
{
    m_gl_tex = new QOpenGLTexture(QOpenGLTexture::Target2D);
}

void RenderLayer::draw(GLuint empty_vao, bool paths_display_active)
{
    QMutexLocker locker(&m_mutex);

    if (m_refresh_gl_texture)
    {
        m_gl_tex->destroy();
        m_gl_tex->setMinMagFilters(QOpenGLTexture::Linear, QOpenGLTexture::Nearest);
        m_gl_tex->setData(m_image, QOpenGLTexture::MipMapGeneration::DontGenerateMipMaps);
        m_refresh_gl_texture = false;
    }

    m_gl->glUseProgram(m_shader_program);

    GLfloat mult = paths_display_active ? 0.6 : 1.0;
    m_gl->glUniform1f(m_mult_loc, mult);

    m_gl->glActiveTexture(GL_TEXTURE0);
    m_gl_tex->bind();
    m_gl->glDisable(GL_DEPTH_TEST);
    m_gl->glDepthMask(GL_FALSE);
    m_gl->glBindVertexArray(empty_vao);
    m_gl->glDrawArrays(GL_TRIANGLES, 0, 3);
    m_gl->glDepthMask(GL_TRUE);
}

void RenderLayer::init_gl(QSurfaceFormat format)
{
    if (!m_gl)
    {
        RENDERER_LOG_ERROR("Attempted to initialize GL without first setting GL functions");
        return;
    }

    auto vertex_shader = load_gl_shader("fullscreen_tri.vert");
    auto fragment_shader = load_gl_shader("final_render.frag");

    m_shader_program = create_shader_program(
        m_gl,
        &vertex_shader,
        &fragment_shader);

    m_mult_loc = m_gl->glGetUniformLocation(m_shader_program, "u_mult");

    m_gl_initialized = true;
}

void RenderLayer::set_gl_functions(QOpenGLFunctions_4_1_Core* functions)
{
    m_gl = functions;
}

void RenderLayer::set_display_transform(const QString& transform)
{
    // todo: redundant with RenderWidget::slot_display_transform_changed.

    QMutexLocker locker(&m_mutex);

    OCIO::DisplayTransformRcPtr transform_ptr = OCIO::DisplayTransform::Create();
    transform_ptr->setInputColorSpaceName(OCIO::ROLE_SCENE_LINEAR);
    transform_ptr->setDisplay(m_ocio_config->getDefaultDisplay());
    transform_ptr->setView(transform.toStdString().c_str());

    OCIO::ConstContextRcPtr context = m_ocio_config->getCurrentContext();
    m_ocio_processor = m_ocio_config->getProcessor(context, transform_ptr, OCIO::TRANSFORM_DIR_FORWARD);

    if (m_image_storage)
    {
        const CanvasProperties& frame_props = m_image_storage->properties();
        for (std::size_t y = 0; y < frame_props.m_tile_count_y; ++y)
        {
            for (std::size_t x = 0; x < frame_props.m_tile_count_x; ++x)
                update_tile_no_lock(x, y);
        }
    }
}

}   // namespace studio
}   // namespace appleseed
