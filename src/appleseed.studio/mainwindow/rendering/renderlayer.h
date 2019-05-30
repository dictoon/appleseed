
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

#pragma once

// appleseed.studio headers.
#include "mainwindow/rendering/renderwidget.h"

// OpenColorIO headers.
#include <OpenColorIO/OpenColorIO.h>
namespace OCIO = OCIO_NAMESPACE;

// Qt headers.
#include <QOpenGLWidget>

// Standard headers.
#include <cstddef>

// Forward declarations.
class QOpenGLFunctions_4_1_Core;
class QOpenGLTexture;
class QWidget;

namespace appleseed {
namespace studio {

//
// A render widget based on QImage.
//

class RenderLayer
  : public RenderWidget
{
    Q_OBJECT

  public:
    // Constructor.
    RenderLayer(
        const std::size_t       width,
        const std::size_t       height,
        OCIO::ConstConfigRcPtr  ocio_config,
        QWidget*                parent = nullptr);

    // Thread-safe.
    void set_display_transform(
        const QString&          transform);

    void draw(
        const GLuint            empty_vao,
        const bool              paths_display_active);

    void init_gl(QSurfaceFormat format);
    void set_gl_functions(QOpenGLFunctions_4_1_Core* functions);

  private:
    QOpenGLFunctions_4_1_Core*  m_gl;
    QOpenGLTexture*             m_gl_tex;
    GLuint                      m_shader_program;
    GLint                       m_mult_loc;
    bool                        m_gl_initialized;
    bool                        m_refresh_gl_texture;
};

}   // namespace studio
}   // namespace appleseed
