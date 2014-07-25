

//
// This source file is part of appleseed.
// Visit http://appleseedhq.net/ for additional information and resources.
//
// This software is released under the MIT license.
//
// Copyright (c) 2010-2013 Francois Beaune, Jupiter Jazz Limited
// Copyright (c) 2014 Francois Beaune, The appleseedhq Organization
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

#ifndef APPLESEED_RENDERER_MODELING_INPUT_EXPRESSION_H
#define APPLESEED_RENDERER_MODELING_INPUT_EXPRESSION_H

// appleseed.renderer headers.
#include "renderer/global/globaltypes.h"

// appleseed.foundation headers.
#include "foundation/image/color.h"
#include "foundation/platform/compiler.h"
#include "foundation/platform/types.h"
#include "foundation/core/concepts/noncopyable.h"

// SeExpr headers
#include "SeExpression.h"

// Forward declarations.
namespace renderer { class ShadingPoint; }

namespace renderer
{

//
// Expression.
//

class DLLSYMBOL Expression
{
  public:
    // Constructor.
    Expression();
    
    // Constructor.
    explicit Expression(const char* expr, bool is_vector = true);

    // Destructor.
    ~Expression();

    // Copy constructor.
    Expression(const Expression& other);

    // Assignment.
    Expression& operator=(const Expression& other);

    void swap(Expression& other);

    void set_expression(const char* expr, bool is_vector = true);

    bool syntax_ok() const;

    foundation::Color3d evaluate(const ShadingPoint& shading_point) const;

  private:
    struct Impl;
    Impl *impl;
};

}       // namespace renderer

#endif  // !APPLESEED_RENDERER_MODELING_INPUT_EXPRESSION_H
