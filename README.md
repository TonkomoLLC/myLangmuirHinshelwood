# myLangmuirHinshelwood

Example of how to compile your own reaction rate expression.

myLangmuirHinshelwood functions exactly like the LangmuirHinshelwood reaction rate in OpenFOAM 12. 

This code will compile libmyLangmuriHinshelwood.so, which can be added to a case and then edited for your own custom reaction rate.

# Rate Equation 

The rate equation for the `irreversibleMyLangmuirHinshelwoodReaction` is given by:

$$
\text{Rate} = \frac{A_0 \, T^{\beta_0} \, e^{-T_{a,0}/T}}{\left[\, a + A_1 \, T^{\beta_1} \, e^{-T_{a,1}/T} \, C_{\text{CH4}}^{m_1} + A_2 \, T^{\beta_2} \, e^{-T_{a,2}/T} \, C_{\text{O2}}^{m_2} \,\right]^{m_0}}
$$

Where:

- $A_i$, $T_{a,i}$, $\beta_i$, and $m_i$ are the pre-exponential factors, activation temperatures, temperature exponents, and reaction orders for each term.
- $a$ is the adsorption constant.
- $C_{\text{CH4}}$ and $C_{\text{O2}}$ are the concentrations of methane and oxygen.

---

## Parameter Table

| Parameter | Description | Units | Example Values |
| :-- | :-- | :-- | :-- |
| $A$ | Pre-exponential factors | [m³/kmol/s] | $(5.2 \times 10^{16}, 1.0 \times 10^{10}, 2.0 \times 10^8)$ |
| $T_a$ | Activation temperatures | [K] | $(14906, 10000, 8000)$ |
| $\beta$ | Temperature exponents | - | $(0, 0, 0)$ |
| $m$ | Reaction orders | - | $(2, 1, 1)$ |
| $a$ | Adsorption constant | - | $1.0$ |

---

## Notes

- The denominator models competitive adsorption effects.
- The `reactants` list in the reaction dictionary specifies the species order (e.g., `reactants (CH4 O2)`).
- The arrays for `A`, `Ta`, `beta`, and `m` must contain exactly **three values** corresponding to the three terms in the rate equation.



# How to modify for your own reaction rate

- read the proper dictionary entries

- modify the member functions

   - `operator()`: Computes the reaction rate itself.
   - `ddc`: Computes derivatives with respect to species concentrations.
   - `ddT`: Specifically handles temperature sensitivity.



## ddc

the `ddc` function computes the **derivative of the reaction rate with respect to concentration** for each species involved in the reaction.

Specifically, the function signature:

```cpp
inline void ddc(
    const scalar p,
    const scalar T,
    const scalarField&amp; c,
    const label li,
    scalarField&amp; ddc
) const;
```

is designed to fill the `ddc` array with the partial derivatives of the reaction rate with respect to each species' concentration (i.e., $\frac{\partial \text{rate}}{\partial c_i}$).

Within the implementation, you can see:

- `ddc` is initialized to zero for all species.
- The relevant derivatives are then assigned to the indices corresponding to the reactant species:

```cpp
ddc[r_[^0]] = dkdb*dbdc1;
ddc[r_[^1]] = dkdb*dbdc2;
```

where `dkdb*dbdc1` and `dkdb*dbdc2` represent the analytical derivatives of the rate with respect to the concentrations of the first and second reactants, respectively[^2].

This matches the typical OpenFOAM convention for reaction rate classes, where `ddc` provides the vector of partial derivatives of the rate with respect to the concentration field. The function `hasDdc()` returning `true` further confirms that this reaction rate supports and provides these derivatives


## ddT

The `ddT` function computes the **derivative of the reaction rate with respect to temperature** (\$ \frac{\partial rate}{\partial T} \$), not `ddt` (which typically refers to time derivatives in PDE contexts). 

---

### Key Implementation Details

1. **Function Signature**

```cpp
inline scalar ddT(
    const scalar p,
    const scalar T,  // Temperature
    const scalarField&amp; c,
    const label li
) const;
```

Explicitly takes temperature `T` as an input and returns a `scalar` derivative value.
2. **Mathematical Formulation**
The derivative calculation includes terms like:

```cpp
const scalar dk0dT = k0/T*(beta_[^0] + Ta_[^0]/T);
const scalar dk1dT = k1/T*(beta_[^1] + Ta_[^1]/T);
const scalar dk2dT = k2/T*(beta_[^2] + Ta_[^2]/T);
```

These terms derive from the **Arrhenius-type temperature dependence** (\$ A T^\beta e^{-T_a/T} \$)[^2].

3. **Return Value**

The final expression:

```cpp
return (dk0dT - k0*m_[^0]*dbdT/b)/pow(b, m_[^0]);
```

Represents $\frac{\partial}{\partial T} \left( \frac{k_0}{b^{m_0}} \right)$, where $b$ depends on temperature through $ k_1(T)$ and $ k_2(T)$.


