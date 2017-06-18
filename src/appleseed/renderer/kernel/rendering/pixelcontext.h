
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

#ifndef APPLESEED_RENDERER_KERNEL_RENDERING_PIXELCONTEXT_H
#define APPLESEED_RENDERER_KERNEL_RENDERING_PIXELCONTEXT_H

// appleseed.foundation headers.
#include "foundation/math/vector.h"

// Standard headers.
#include <cstddef>

namespace renderer
{

//
// This class identifies the pixel and sample currently being rendered throughout the
// Tile Renderer -> Pixel Renderer -> Sample Renderer -> Surface Shader -> Lighting Engine chain.
//

class PixelContext
{
  public:
    // Constructor.
    PixelContext(
        const foundation::Vector2i& pixel_coords,
        const size_t                sample_id,
        const foundation::Vector2d& sample_position);

    // Return the coordinates of the pixel being rendered.
    const foundation::Vector2i& get_pixel_coords() const;

    // Return the sample number within the pixel being rendered.
    size_t get_sample_id() const;

    // Return the coordinates (in NDC) of the sample being rendered.
    const foundation::Vector2d& get_sample_position() const;

  private:
    const foundation::Vector2i      m_pixel_coords;
    const size_t                    m_sample_id;
    const foundation::Vector2d      m_sample_position;
};


//
// PixelContext class implementation.
//

inline PixelContext::PixelContext(
    const foundation::Vector2i&     pixel_coords,
    const size_t                    sample_id,
    const foundation::Vector2d&     sample_position)
  : m_pixel_coords(pixel_coords)
  , m_sample_id(sample_id)
  , m_sample_position(sample_position)
{
}

inline const foundation::Vector2i& PixelContext::get_pixel_coords() const
{
    return m_pixel_coords;
}

inline size_t PixelContext::get_sample_id() const
{
    return m_sample_id;
}

inline const foundation::Vector2d& PixelContext::get_sample_position() const
{
    return m_sample_position;
}

}       // namespace renderer

#endif  // !APPLESEED_RENDERER_KERNEL_RENDERING_PIXELCONTEXT_H
