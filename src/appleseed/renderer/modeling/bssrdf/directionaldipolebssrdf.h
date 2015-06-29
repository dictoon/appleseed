
//
// This source file is part of appleseed.
// Visit http://appleseedhq.net/ for additional information and resources.
//
// This software is released under the MIT license.
//
// Copyright (c) 2015 Esteban Tovagliari, The appleseedhq Organization
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

#ifndef APPLESEED_RENDERER_MODELING_BSSRDF_DIRECTIONALDIPOLEBSSRDF_H
#define APPLESEED_RENDERER_MODELING_BSSRDF_DIRECTIONALDIPOLEBSSRDF_H

// appleseed.renderer headers.
#include "renderer/global/globaltypes.h"
#include "renderer/modeling/bssrdf/bssrdf.h"
#include "renderer/modeling/bssrdf/ibssrdffactory.h"

// appleseed.foundation headers.
#include "foundation/math/basis.h"
#include "foundation/math/vector.h"

// appleseed.main headers.
#include "main/dllsymbol.h"

// Standard headers.
#include <cstddef>

// Forward declarations.
namespace foundation    { class IAbortSwitch; }
namespace renderer      { class Assembly; }
namespace renderer      { class InputEvaluator; }
namespace renderer      { class ParamArray; }
namespace renderer      { class Project; }
namespace renderer      { class ShadingContext; }
namespace renderer      { class ShadingPoint; }

namespace renderer
{

//
// Directional dipole BSSRDF input values.
//

APPLESEED_DECLARE_INPUT_VALUES(DirectionalDipoleBSSRDFInputValues)
{
    Spectrum    m_reflectance;
    Spectrum    m_mean_free_path;
    double      m_mean_free_path_multiplier;
    double      m_anisotropy;
    double      m_from_ior;
    double      m_to_ior;

    // Precomputed values.
    Spectrum    m_cdf;
    double      m_max_mean_free_path;

    // Convenience accessors for use after reflectance and
    // mean free path are converted to absorption and scattering coefficients.

    const Spectrum& sigma_a() const
    {
        return m_reflectance;
    }

    Spectrum& sigma_a()
    {
        return m_reflectance;
    }

    const Spectrum& sigma_s_prime() const
    {
        return m_mean_free_path;
    }

    Spectrum& sigma_s_prime()
    {
        return m_mean_free_path;
    }
};


//
// Directional dipole BSSRDF factory.
//

class APPLESEED_DLLSYMBOL DirectionalDipoleBSSRDFFactory
  : public IBSSRDFFactory
{
  public:
    // Return a string identifying this BSDF model.
    virtual const char* get_model() const APPLESEED_OVERRIDE;

    // Return metadata for this BSDF model.
    virtual foundation::Dictionary get_model_metadata() const APPLESEED_OVERRIDE;

    // Return metadata for the inputs of this BSDF model.
    virtual foundation::DictionaryArray get_input_metadata() const APPLESEED_OVERRIDE;

    // Create a new BSSRDF instance.
    virtual foundation::auto_release_ptr<BSSRDF> create(
        const char*         name,
        const ParamArray&   params) const APPLESEED_OVERRIDE;
};

}       // namespace renderer

#endif  // !APPLESEED_RENDERER_MODELING_BSSRDF_DIRECTIONALDIPOLEBSSRDF_H
