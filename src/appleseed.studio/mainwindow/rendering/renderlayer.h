
//
// This source file is part of appleseed.
// Visit https://appleseedhq.net/ for additional information and resources.
//
// This software is released under the MIT license.
//
// Copyright (c) 2019-2020 Gray Olson, The appleseedhq Organization
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

// appleseed.qtcommon headers.
#include "widgets/renderwidget.h"

// OpenColorIO headers.
#include <OpenColorIO/OpenColorIO.h>
namespace OCIO = OCIO_NAMESPACE;

// Qt headers.
#include <QObject>
#include <QOpenGLTexture>
#include <QOpenGLWidget>

// Standard headers.
#include <cstddef>
#include <memory>

// Forward declarations.
class QOpenGLFunctions_4_1_Core;
class QWidget;

namespace appleseed {
namespace studio {

//
// Turn a render widget into a render layer.
//

class RenderLayer
  : public qtcommon::RenderWidget
{
    Q_OBJECT

  public:
    // Constructor.
    RenderLayer(
        const std::size_t       width,
        const std::size_t       height,
        OCIO::ConstConfigRcPtr  ocio_config,
        QWidget*                parent = nullptr);

    void set_gl_functions(QOpenGLFunctions_4_1_Core* functions);
    void init_gl();

    void draw(
        const GLuint            empty_vao,
        const bool              paths_display_active);

  private:
    QOpenGLFunctions_4_1_Core*          m_gl = nullptr;
    std::unique_ptr<QOpenGLTexture>     m_gl_tex;
    GLuint                              m_shader_program = 0;
    GLint                               m_mult_location = 0;
    bool                                m_refresh_gl_texture = true;
};

}   // namespace studio
}   // namespace appleseed
