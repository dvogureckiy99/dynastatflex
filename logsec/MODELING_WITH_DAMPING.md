
- 2.2 Model Reduction
		- thick rod dynamics.
			- толстый стержень
		- _geometrically exact beam theory_
		- _ﬁnite strain beam formulation_
		- strain
			- напряжение ср
			  нагрузка ж
			  натяжение ср
			  напряженность ж
		- The strong form of the Simo–Reissner beam theory
			- 6 PDEs
			-
			- $\large \bm{f}^{\text{internal}}$ and $\large \bm{m}^{\text{internal}}$ results from internal stresses acting on the beam cross-section area at point r of the beam center line.
			- $\large \bm{f}^{\text{ext}}$ and $\large \bm{m}^{\text{ext}}$ externally imposed **distributed** forces [N/m],[N]
			- $\large \bm{f}^{\text{inertia}}$ and $\large \bm{m}^{\text{inertia}}$ inertia effects
			- the detailed constitutive equations that relate $\large \bm{f}^{\text{internal}}$ and $\large \bm{m}^{\text{internal}}$ to the ﬁrst Piola–Kirchhoff stress tensor require an introduction into 3D continuum mechanics and is omitted in this work.
				- Piola–Kirchhoff stress tensor
				- deﬁne the expressions only after reduction to the speciﬁed special cases
			- We further assume a hyperelastic constitutive relation between these kinetic and kinematic quantities
	- Assumption: Vanishing Shear Strains (Kirchhoff–Love Beam Theory)
		- **Neglecting [shear deformations](((659582c4-d3ac-4557-8cd6-b984bfbdac39)))**
		- **assuming that the crosssection is always perpendicular to the center line of the beam**
		- $\large \bm{f}^{\text{internal}}$ split up to into a $\large {}^{\parallel}(\bm{f}^{\text{internal}})_l$ and $\large {}^{\bot}(\bm{f}^{\text{internal}})_l$
		- the moment balance Eq. 1b reduces to the projection onto the center line tangential base vector $\hat{\bm{g}}_1$
		- Kirchhoff–Love beam equations
			- ![image.png](../assets/image_1704312317392_0.png){:height 182, :width 537}
			- Why $\large \bm r\cdot f$ neglected
			- Eq. 2b is now a scalar expression
	- Assumption: Initially Straight and Isotropic
		- **initially straight beam with an isotropic cross-section**
		- isotropic cross-section
			- (of an object or substance) having a physical property which has the same value when measured in different directions
		- ![image.png](../assets/image_1704312804410_0.png){:height 200, :width 532}
		- $$\large M_{damper}=B_z\cdot (\bm{\kappa})_{t}=\mu \mathsf{w} (\bm{\kappa}) \cdot (\bm{\kappa})_{t}$$
		- $$\large (M_{damper})_l=(B_z)_l\cdot (\bm{\kappa})_{t}+B_z\cdot (\bm{\kappa})_{tl}=\mu \mathsf{w}( (\bm{\kappa})_{l}(\bm{\kappa})_{t} + (\bm{\kappa}) (\bm{\kappa})_{tl} )$$
		- $$\large M_{spring}=EI(\bm{\kappa})$$
		- $$ \begin{align*} &{}^{\parallel}(\bm{f}^{\text{internal}})_l+\left( \frac{(\bm{r})_l}{\parallel(\bm{r})_l\parallel_2^2}\times \left( \underbrace{EI(\bm{\kappa})_l+\mu \mathsf{w}( (\bm{\kappa})_{ll}(\bm{\kappa})_{lt} + (\bm{\kappa})_{l} (\bm{\kappa})_{ltl} )}_{ (\bm{m}^{\text{internal}})_l }+\bm{m}^{ext}+\bm{m}^{inertia} \right) \right)_l- &&(3a)\\ &- \bm{f}^{ext} \underbrace{- \rho A (\bm{r})_{tt}}_{\bm{f}^{inertia}} = \bm{0} && \end{align*}$$
		- ![image.png](../assets/image_1705421701886_0.png)
		- Young modulus E, inertia I, density ρ, cross-section area A
		- this all wrong!
		  background-color:: red
			- density there is **area density**, i.e. $\rho=\rho_V\cdot dl$, where  $\rho_V$ — **volumetric mass density**, $dl$ — elementary beam length.
			  background-color:: red
			- $dl=\Delta l$ as we're making FEM formulation, where $\Delta l$ is length of finite elements.
			  background-color:: red
			  and $\rho=\rho_V\cdot \Delta l$
			- $\large 2\rho I \omega_t$            this is            $\Large \frac{kg}{m^2}m^4\frac{1}{s}=\frac{kg*m^2}{s}\,\, [M]$
			  background-color:: red
			- better amd more correct would be there express $\rho A$ as $m_{d l}=\rho_V\cdot d l\cdot\mathsf{w}\cdot h$ or after FEM formulation $m_{\Delta l}=\rho_V\cdot \Delta l\cdot\mathsf{w}\cdot h$
			  background-color:: red
			- $\mathsf{w}$ - width of side that's bending
			  background-color:: red
			- $h$ - height of side that's not influencing on bending
			  background-color:: red
			-
		- curvature vector
		  collapsed:: true
			- $\Large \kappa=\frac{y^{\prime \prime} x^{\prime}-x^{\prime \prime} y^{\prime}}{\left(x^{\prime2}+y^{\prime2}\right)^{\frac{2}{3}}}$
			  Let **γ**(*t*) = (*x*(*t*), *y*(*t*)) be a proper [parametric representation](https://en.wikipedia.org/wiki/Parametric_representation) of a twice differentiable plane curve.
			  where primes refer to derivatives with respect to *t*.
			- These can be expressed in a coordinate-free way as
			  ![image.png](../assets/image_1704402253032_0.png)
			- For a curve defined by an implicit equation F(x, y) = 0 with partial derivatives denoted Fx , Fy , Fxx , Fxy , Fyy , the curvature is given by
			  ![image.png](../assets/image_1704402450530_0.png)
			- {{embed ((65971a9b-9762-49a6-900a-f8c0d616bf01))}}
			- Normal vector or curvature vector
				- ![image.png](../assets/image_1704312822866_0.png)
				- [from article](((659722e9-22dd-44de-99fc-cb8c87c68479)))
				- A curve *[normal vector](https://en.wikipedia.org/wiki/Normal_vector)*, sometimes called the **curvature vector**, indicates the deviance of the curve from being a straight line.
				  It is defined as
				- ![image.png](../assets/image_1704403056646_0.png)
				- Its normalized form, the unit normal vector, is the second Frenet vector **e**2(*t*) and is defined as
				- ![image.png](../assets/image_1704403076582_0.png)
				- ![image.png](../assets/image_1704403246155_0.png)
				- ![image.png](../assets/image_1704403297163_0.png)
				- The unit tangent vector determines the orientation of the curve, or the 
				  forward direction, corresponding to the increasing values of the 
				  parameter.
		- ![image.png](../assets/image_1704312884180_0.png)
	- Assumption: Torsion-Free
		- **no [torsional effects](((659582f4-2624-4e09-b75d-b83f4ce9e0a9)))**
		- only the perpendicular component of the external moment affects the force balance equation
			- ![image.png](../assets/image_1704313642060_0.png)
		- ![image.png](../assets/image_1704313725218_0.png)
	- Assumption: Inextensible Beam
		- ![image.png](../assets/image_1704314541098_0.png)
		- and thus ${}^{\parallel}(\bm{f}^{\text{internal}})_l$ vanishes.
		- ![image.png](../assets/image_1704314589615_0.png)
		- $$\large \begin{align*} &\left( (\bm{r})_l\times \left( EI(\bm{\kappa})_l+\mu \mathsf{w}( (\bm{\kappa})_{l}(\bm{\kappa})_{t} + (\bm{\kappa}) (\bm{\kappa})_{tl} )+\bm{m}^{ext} \right) \right)_l - \bm{f}^{ext} - \rho A (\bm{r})_{tt} = \bm{0} &&(8b) \end{align*}$$
		- it is in general difﬁcult to ﬁnd such a set of variables that fulﬁll the inextensibility constraint Eq. 8a by construction
			- A common practice to enforce the equality constraint Eq. 8a on the simulation result in a weak sense
				- in an integral form instead of point-wise, is by means of extending the weak form of the model Eq. 8b with a Lagrange multiplier potential
		- The following assumption of pure planar bending, however, does again permit a parametrization that fulﬁlls this constraint directly in the strong sense, i.e., for every point along the beam.
	- Assumption: Pure Planar Bending
		- ((6597d0d0-00bc-46ce-b5a8-f17b1c88ba87))
		- switch to a component-wise notation in Cartesian coordinates
		- ![image.png](../assets/image_1704315328205_0.png)
		- using this:
			- ![image.png](../assets/image_1705421921439_0.png)
		- with this assumptions (4) becomes
			- $\bm{\kappa}=\begin{bmatrix}0& 0& \varphi_l \end{bmatrix}^T$
			- ![image.png](../assets/image_1704316227771_0.png)
				- ![image.png](../assets/image_1704316241037_0.png)
				  id:: 659fb250-5912-4b7a-ba3e-6c67f05645c4
		- ![image.png](../assets/image_1704316386525_0.png)
		- and (8b) becoming (let $\mu \mathsf{w}$ be $\mu^*$)
		- $$ \begin{align*} &\begin{bmatrix} sin(\varphi)(EI(\varphi)_{ll}+\mu^* ( (\varphi)_{ll}(\varphi)_{lt} + (\varphi)_{l} (\varphi)_{ltl} )+{}^{\bot}m^{ext}_z) \\ -cos(\varphi)(EI(\varphi)_{ll}+\mu^* ( (\varphi)_{ll}(\varphi)_{lt} + (\varphi)_{l} (\varphi)_{ltl} )+{}^{\bot}m^{ext}_z) \end{bmatrix}_l + \begin{bmatrix} f^{ext}_x\\f^{ext}_y \end{bmatrix}  - m_{d l} \begin{bmatrix} x\\ y \end{bmatrix}_{tt} && (12)\end{align*}$$
		-
	- 2.2.1 Static Beam Model Expressed in the Curve Tangent Angle
		- ![image.png](../assets/image_1704323361120_0.png)
		- ![image.png](../assets/image_1704323421121_0.png)
		- ![image.png](../assets/image_1704323674098_0.png)
		- multiply 14 by 15 give, i.e. rotation around z-axis by $\varphi$ in clockwise direction to the inertial coordinate system
		- ![image.png](../assets/image_1704324056023_0.png)
		- In case of an absent external force $\large \bm{f}^{\text{ext}}$, the static beam model Eq. 16 even admits a simple analytic solution. If a nontrivial curvature $(\varphi)_l \neq 0$ is assumed, Eq. 16 reduces to
		- ![image.png](../assets/image_1704324575684_0.png){:height 76, :width 360}
		- which can be integrated twice and yields a unique solution if boundary conditions are applied.
	- 2.2.2 Dynamic Beam Model Expressed in the Curve Tangent Angle
		- ((659d8944-78cc-42ef-8d40-62116747cd99))
		- Cartesian xy acceleration terms in Eq. 12 remain to be expressed in terms of the curve tangent angle φ.
		- Assuming no buckling of the object, x and y have continuous derivatives, and thus, **Schwarz’s theorem** allows changing the order of the derivations.
			- Schwarz's theorem let us
				- $$\frac{\partial \partial^2 \bm{r}(l,t)}{\partial l\partial t^2}=\frac{\partial^2 \partial \bm{r}(l,t)}{\partial t^2 \partial l}$$
			- ((659ee7be-7fd0-4f87-8df6-d321b2a3a782)) Why this is need?
				- ![image.png](../assets/image_1704913222851_0.png)
			- $$\large \begin{align*} &\begin{bmatrix} sin(\varphi)\varphi^* \\ -cos(\varphi)\varphi^* \end{bmatrix}_{ll} + \begin{bmatrix} f^{ext}_x\\f^{ext}_y \end{bmatrix}_l  - m_{d l} \begin{bmatrix} cos(\varphi)\\ sin(\varphi) \end{bmatrix}_{tt}=\bm{0} && (18a)\end{align*}$$
			  
			  $$\large \begin{align*} &\varphi^*=EI(\varphi)_{ll}+\mu^* ( (\varphi)_{ll}(\varphi)_{lt} + (\varphi)_{l} (\varphi)_{ltl} )+{}^{\bot}m^{ext}_z &&(18b) \end{align*}$$
		- PDE system entirely expressed in terms of the curve tangent angle φ.
		- As for the static case (Eq. 16), the geometric identities fulﬁll the inextensibility constraint Eq. 9b by construction; thus, no special consideration is necessary.
		- expanding all partial derivatives and grouping the trigonometric terms, Eq. 18a yields
			- $$\large \begin{bmatrix} cos(\varphi)\varphi_l\varphi^*+sin(\varphi)\varphi^*_{l} \\ 
			  sin(\varphi)\varphi_l\varphi^*-cos(\varphi)\varphi^*_{l} \end{bmatrix}_l$$
			- $$\large \begin{bmatrix} (-sin(\varphi)\varphi^2_l+cos(\varphi)\varphi_{ll})\varphi^*+cos(\varphi)\varphi_l\varphi^*_{l}+cos(\varphi)\varphi_l\varphi^*_{l}+sin(\varphi)\varphi^*_{ll} \\ 
			  (cos(\varphi)\varphi^2_l+sin(\varphi)\varphi_{ll})\varphi^*+sin(\varphi)\varphi_l\varphi^*_{l}+sin(\varphi)\varphi_l\varphi^*_{l}-cos(\varphi)\varphi^*_{ll} \end{bmatrix}$$
			- $$\large \begin{bmatrix} cos(\varphi)[\varphi_{ll}\varphi^*+2\varphi_l\varphi^*_{l}]-sin(\varphi)[\varphi^2_l\varphi^*-\varphi^*_{ll}] \\ sin(\varphi)[\varphi_{ll}\varphi^*+2\varphi_l\varphi^*_{l}]+cos(\varphi)[\varphi^2_l\varphi^*-\varphi^*_{ll}] \end{bmatrix}$$
			- Similar to the static case, premultiplying the entire system with the rotation matrix $\bm{R}_z(\varphi)$ from Eq. 15 again extracts the components parallel and perpendicular to the beam center line.
			- $$\large \begin{bmatrix} cos(\varphi) & sin(\varphi) \\ -sin(\varphi) & cos(\varphi) \end{bmatrix}\cdot\begin{bmatrix} cos(\varphi)a - sin(\varphi)b \\ sin(\varphi)a+cos(\varphi)b \end{bmatrix}=\begin{bmatrix}a\\b\end{bmatrix}$$
			- $$\large \begin{bmatrix} [\varphi_{ll}\varphi^*+2\varphi_l\varphi^*_{l}] \\ [\varphi^2_l\varphi^*-\varphi^*_{ll}] \end{bmatrix}$$
			- $- \mu \bm{R}_z(\varphi) \begin{bmatrix} cos(\varphi)\\sin(\varphi) \end{bmatrix}_t=- \mu \bm{R}_z(\varphi) \begin{bmatrix} -sin(\varphi)\varphi_t\\cos(\varphi)\varphi_t \end{bmatrix}= -\mu \begin{bmatrix} -cos(\varphi)sin(\varphi)\varphi_t+sin(\varphi)cos(\varphi)\varphi_t \\ sin^2(\varphi)\varphi_t+cos^2(\varphi)\varphi_t\end{bmatrix}=-\mu \begin{bmatrix} 0 \\ \varphi_t\end{bmatrix}$
			- all together from (18)
			- $$\large \begin{align*} &\begin{bmatrix} \varphi_{ll}\varphi^*+2\varphi_l\varphi^*_{l}+ m_{dl} \varphi^2_t+({}^{\parallel}f^{ext})_l \\ \varphi^2_l\varphi^*-\varphi^*_{ll}- m_{dl}(\varphi)_{tt}+({}^{\bot}f^{ext})_l \end{bmatrix}=\bm{0} && (19) \end{align*} $$
			- The only acceleration term $(φ)_{tt}$ , however, appears solely in the perpendicular direction
			- $$\large \begin{align*}& && (19a) \end{align*}$$
			- $$\large \begin{align*}&(\varphi^*)_{ll}=EI(\varphi)_{llll}+({}^{\bot}m^{ext}_z)_{ll}+\mu^*( (\varphi)_{ll}(\varphi)_{lt} + (\varphi)_{l} (\varphi)_{ltl} )_{ll}= &&(19a) \\ &=EI(\varphi)_{llll}+({}^{\bot}m^{ext}_z)_{ll}+ M^{d} && \end{align*}$$
			- derivative lalculated with sympy package
				- $$M^d=\mu^* \left(\frac{\partial}{\partial l} \varphi{\left(l,t \right)} \frac{\partial^{5}}{\partial t\partial l^{4}} \varphi{\left(l,t \right)} + 3 \frac{\partial^{2}}{\partial l^{2}} \varphi{\left(l,t \right)} \frac{\partial^{4}}{\partial t\partial l^{3}} \varphi{\left(l,t \right)}+\right.\\ \left.+ 3 \frac{\partial^{3}}{\partial l^{3}} \varphi{\left(l,t \right)} \frac{\partial^{3}}{\partial t\partial l^{2}} \varphi{\left(l,t \right)} + \frac{\partial^{4}}{\partial l^{4}} \varphi{\left(l,t \right)} \frac{\partial^{2}}{\partial t\partial l} \varphi{\left(l,t \right)}\right)$$
					- $$\Large M^d=\mu^*(\varphi_l\varphi_{llllt}+3\varphi_{ll}\varphi_{lllt}+3\varphi_{lll}\varphi_{llt}+\varphi_{llll}\varphi_{lt})$$
			- if $\varphi^M=\begin{bmatrix} \varphi_l & \sqrt{3}\varphi_{ll} & \sqrt{3}\varphi_{lll}&\varphi_{llll} \end{bmatrix}^T$, then
				- $$M^d=\mu^* (\varphi^M)^T \times \text{flip}((\varphi^M)_t)$$
				  $$(\varphi^*)_{ll}=EI(\varphi)_{llll}+({}^{\bot}m^{ext}_z)_{ll}+ M^{d}$$
			- if $\varphi^{M1}=\begin{bmatrix} \varphi_l &\varphi_{ll} \end{bmatrix}^T$, then
				- $$\varphi^*=EI(\varphi)_{ll}+({}^{\bot}m^{ext}_z)+M^{d1}$$
				  $$M^{d1}=\mu^* (\varphi^{M1})^T \times \text{flip}((\varphi^{M1})_t)$$
			- then all together
				- $$\varphi^2_l[EI(\varphi)_{ll}+({}^{\bot}m^{ext}_z)+M^{d1}]-[EI(\varphi)_{llll}+({}^{\bot}m^{ext}_z)_{ll}+ M^{d}]- m_{dl}(\varphi)_{tt}+({}^{\bot}f^{ext})_l =0$$
				- $$ m_{dl}(\varphi)_{tt}=EI[\varphi^2_l (\varphi)_{ll} - (\varphi)_{llll} ]+ [\varphi^2_l M^{d1} - M^{d}]+\\+({}^{\bot}f^{ext})_l + [\varphi^2_l{}^{\bot}m^{ext}_z-({}^{\bot}m^{ext}_z)_{ll}]$$
				- $$\large m_{dl}(\varphi)_{tt} + M^{d}(\varphi)-\varphi^2_l M^{d1}(\varphi) =EI[\varphi^2_l (\varphi)_{ll} - (\varphi)_{llll} ]+ \\+({}^{\bot}f^{ext})_l + [\varphi^2_l{}^{\bot}m^{ext}_z-({}^{\bot}m^{ext}_z)_{ll}]$$(20a)
				- $$\large (\varphi)_{tt} + c_2 ( M^{d}_{/\mu^*}(\varphi)-\varphi^2_l M^{d1}_{/\mu^*}(\varphi) )=c_1( [\varphi^2_l (\varphi)_{ll} - (\varphi)_{llll} ]+ \\+c_3({}^{\bot}f^{ext})_l + c_3[\varphi^2_l{}^{\bot}m^{ext}_z-({}^{\bot}m^{ext}_z)_{ll}])$$(20b)
				  with $\Large c_1 = \frac{EI}{m_{dl}}=\frac{EI}{\rho A}$ and $\Large c_2=\frac{\mu^*}{\rho A}=\frac{\mu \mathsf{w}}{\rho A}$ and $\Large c_3=\frac{1}{\rho A}$
			- which in the case of no external inputs admits the very concise strong form
				- $$\large \begin{align*} &(\varphi)_{tt} + c_2 ( M^{d}_{/\mu^*}(\varphi)-\varphi^2_l M^{d1}_{/\mu^*}(\varphi) )=c_1( [\varphi^2_l (\varphi)_{ll} - (\varphi)_{llll} ] &&(21) \end{align*}$$
				- $$\large \begin{align*} &(\varphi)_{tt} + c_2 ( (\varphi^{M})^T \times \text{flip}((\varphi^{M})_t)-\varphi^2_l (\varphi^{M1})^T \times \text{flip}((\varphi^{M1})_t))= &&(21b)\\&=c_1( [\varphi^2_l (\varphi)_{ll} - (\varphi)_{llll} ] && \end{align*}$$
				- $$\begin{align*} &(\varphi)_{tt} + c_2 [ (\varphi_l\varphi_{llllt}+3\varphi_{ll}\varphi_{lllt}+3\varphi_{lll}\varphi_{llt}+\varphi_{llll}\varphi_{lt})-&&\\&- \varphi^2_l ( (\varphi)_{ll}(\varphi)_{lt} + (\varphi)_{l} (\varphi)_{ltl} ) ]= &&(21c)\\&=c_1( [\varphi^2_l (\varphi)_{ll} - (\varphi)_{llll} ] && \end{align*}$$
				- $$\Large \begin{align*} &(\varphi)_{tt} + c_2 [ \varphi_l\varphi_{llllt}+3\varphi_{ll}\varphi_{lllt}+3\varphi_{lll}\varphi_{llt}+\varphi_{llll}\varphi_{lt}-&&\\&- (\varphi)^2_l (\varphi)_{ll}(\varphi)_{lt} - (\varphi)^3_{l} (\varphi)_{llt} ) ]= &&(21d)\\&=c_1( [\varphi^2_l (\varphi)_{ll} - (\varphi)_{llll} ] && \end{align*}$$
			- single PDE governing the beam dynamics in a single parameter φ
			- While this reduced model is relevant for PDE controller development, it is not directly applicable for use in simulations. We therefore present in the following section a respective approximation with a system of ordinary differential equations (ODEs), in terms of a FEM formulation.
- FEM FORMULATION
	- ((6596d79f-7c57-4a3c-b240-610564f9accb))
	- we outline the development of a FEM simulation procedure, starting from the development of the weak form of the beam model (Eq. 21), without considering external forces.
	- After transforming the integrodifferential weak form into a system of nonlinear ODEs of thee second order in time via a Bubnov–Galerkin approximation, a ﬁnite element discretization leads to a simulation procedure
	- **Weak Form of Large Deformation in the Curve Tangent Angle**
		- $$\Large \begin{align*} &\frac{1}{c_1}\int_0^L(\varphi)_{tt}\delta \varphi dl + \frac{c_2}{c_1} \left[ \int_0^L \varphi_l\varphi_{llllt}\delta \varphi dl+3\int_0^L\varphi_{ll}\varphi_{lllt}\delta \varphi dl+\right.&&\\&+3\int_0^L\varphi_{lll}\varphi_{llt}\delta \varphi dl+\int_0^L\varphi_{llll}\varphi_{lt}\delta \varphi dl- \int_0^L(\varphi)^2_l (\varphi)_{ll}(\varphi)_{lt}\delta \varphi dl -&&(22)\\&\left.- \int_0^L(\varphi)^3_{l} (\varphi)_{llt}\delta \varphi dl \right] =\int_0^L\varphi^2_l (\varphi)_{ll}\delta \varphi dl - \int_0^L(\varphi)_{llll}\delta \varphi dl && \end{align*}$$
		- this is the final dynamic equation in a weak form that builds the basis for the following FEM formulation.
		- ![image.png](../assets/image_1705925270164_0.png)
	- **Bubnov–Galerkin Approximation**
		- $$\Large \begin{align*}& \varphi(l,t)\approx\varphi^h(l,t)=\sum_{e_1=1}^N\bm{\xi}_{e_1}^a=\sum_{e_1=1}^N(\bm{\psi}_{e_1}(l))^T\bm{a}_{e_1}^{\varphi}(t), &&(25a)\\&\delta \varphi(l,t)\approx\delta\varphi^h(l,t)=\sum_{e_2=1}^N\bm{\xi}_{e_2}^b=\sum_{e_2=1}^N(\bm{\psi}_{e_2}(l))^T\bm{b}_{e_2}^{\varphi}(t), &&(25b) \end{align*}$$
		- $$\large \begin{align*} &\left\{ \begin{align*} 1 \\ 2 \end{align*} \right. ,&& (19a) \end{align*}$$
		- As the coefficients bφ j from the variation Eq. 25b are arbitrary, the weak formulation result in the system of n equations by substituting consecuently $\{b_1=1,b_{j\neq 1}=0;\,b_2=1,b_{j\neq 2}=0\}$ and so on
		  background-color:: green
		-
		- using the same set of n weighted orthogonal spatial basis functions ψ1...n ∈ H2 , together with n time-dependent scaling coefficients aφ i for then approximation of the curve tangent angle and bφ j for the test function.
		-
		- Reorganizing the terms and using a matrix representation ﬁnally leads to
			-
			- ![image.png](../assets/image_1704732642916_0.png)
			- ![image.png](../assets/image_1704733346726_0.png)
			- $$\Large \begin{align*}&\bm{M}(\bm{a}^{\varphi})_{tt}=-c_1\left[\bm{F}\bm{a}^{\varphi}+\frac{1}{3}\left(\bm{f}^3(\bm{a}^{\varphi})-\bm{F}^3(\bm{a}^{\varphi}) \right)\right], &&(28a) \end{align*}$$
			- $$\Large \begin{align*}&f_j^3(\bm{a}^{\varphi})=\int_0^L\left((\bm{a}^{\varphi})^T(\bm{\psi})_l \right)^3(\psi_j)_l dl=\int_0^L\left(\sum_i a_i^{\varphi}(\psi_i)_l \right)^3(\psi_j)_l dl, &&(28d) \end{align*}$$
			- $$\Large \begin{align*}&F_j^3(\bm{a}^{\varphi})=\left((\bm{a}^{\varphi})^T(\bm{\psi})_l \right)^3\psi_j \Bigg|_0^L, &&(28e) \end{align*}$$
			-
			- ![image.png](../assets/image_1704733307972_0.png)
	- **Finite Element Discretization**
		- $$\large \begin{align*}& \begin{bmatrix} \psi_1(\lambda) \\ \psi_2(\lambda) \\ \psi_3(\lambda)\\ \psi_4(\lambda)\\ \psi_5(\lambda)\\ \psi_6(\lambda) \end{bmatrix} =\begin{bmatrix}
		    -6 & 15 & -10 & 0 & 0 & 1\\
		    -3 & 8 & -6 & 0 & 1 & 0\\
		    -1/2 & 3/2 & -3/2 & 1/2 & 0 & 0\\
		    6 & -15 & 10 & 0 & 0 & 0\\
		    -3 & 7 & -4 & 0 & 0 & 0\\
		    1/2 & -1 & 1/2 & 0 & 0 & 0\\
		  \end{bmatrix} \begin{bmatrix} \lambda^5\\\lambda^4\\\lambda^3\\\lambda^2\\\lambda^1\\\lambda^0 \end{bmatrix},&\lambda\in[0,1].  && (29) \end{align*}$$
		- $$\large \begin{align*}& \bm{\xi}_e \large = \underbrace{ \begin{bmatrix} \psi_1^e(\lambda)&\psi_2^e(\lambda)&\psi_3^e(\lambda)&\psi_4^e(\lambda)&\psi_5^e(\lambda)&\psi_6^e(\lambda) \end{bmatrix}^T}_{\Large \bm{\psi}^T_e(\lambda)} \underbrace{ \begin{bmatrix} (\varphi)(\mathcal{N}_{e-1},t) \\(\varphi)_l(\mathcal{N}_{e-1},t)\\(\varphi)_{ll}(\mathcal{N}_{e-1},t)\\(\varphi)(\mathcal{N}_{e},t) \\(\varphi)_l(\mathcal{N}_{e},t)\\(\varphi)_{ll}(\mathcal{N}_{e},t) \end{bmatrix}}_{\Large \bm{a}_e^{\varphi}(t)}, && (30) \\ &e\in[1,N]\text{ — element's numbers},\,\,\, \mathcal{N}_{0\dots N}\text{ — node's locations}&& \end{align*}$$
		- $$\large \begin{align*} &\psi_i^e(l)=\left\{ \begin{align*} &\psi_i(\lambda) &&\text{,if}&&l\leq\mathcal{N}_e&&and&&l\geq\mathcal{N}_{e-1},&&\lambda=\frac{l-\mathcal{N}_{e-1}}{\Delta l} \\ &0&&\text{,else}&& && && &&\end{align*} \right. ,i\in[1,6] && (30b) \end{align*}$$
		- $$\large \begin{align*} &\psi_m(l)=\left\{ \begin{align*} &\psi_{(m-1)\%6+1}(\lambda) &&\text{,if}&&l\leq\mathcal{N}_{(m-1)//6+1}&&and&&l\geq\mathcal{N}_{(m-1)//6},&& \\ &0&&\text{,else}&& && && &&\end{align*} \right. , && (30c)\\&\lambda=\frac{l-\mathcal{N}_{(m-1)//6}}{\Delta l}, m\in[1,6N] && \end{align*}$$
		- $$\Large \begin{align*} &a_{e,i}^{\varphi}(t)=\left\{ \begin{align*} &(\varphi)(\mathcal{N}_{e-1},t) &&\text{,if}&&i=1\\           &(\varphi)_l(\mathcal{N}_{e-1},t) &&\text{,if}&&i=2\\               &(\varphi)_{ll}(\mathcal{N}_{e-1},t) &&\text{,if}&&i=3\\          &(\varphi)(\mathcal{N}_{e},t) &&\text{,if}&&i=4\\        &(\varphi)_l(\mathcal{N}_{e},t) &&\text{,if}&&i=5\\            &(\varphi)_{ll}(\mathcal{N}_{e},t) &&\text{,if}&&i=6 \end{align*} \right. 
		          ,i\in[1,6] && (30d) \end{align*}$$
		- $$\Large \begin{align*} &a_{m}^{\varphi}(t)=\left\{ \begin{align*} &(\varphi)(\mathcal{N}_{(m-1)//6},t) &&\text{,if}&&m\%6=1\\           &(\varphi)_l(\mathcal{N}_{(m-1)//6},t) &&\text{,if}&&m\%6=2\\               &(\varphi)_{ll}(\mathcal{N}_{(m-1)//6},t) &&\text{,if}&&m\%6=3\\          &(\varphi)(\mathcal{N}_{(m-1)//6+1},t) &&\text{,if}&&m\%6=4\\        &(\varphi)_l(\mathcal{N}_{(m-1)//6+1},t) &&\text{,if}&&m\%6=5\\            &(\varphi)_{ll}(\mathcal{N}_{(m-1)//6+1},t) &&\text{,if}&&m\%6=0 \end{align*} \right.         ,&&m\in[1,6N] && (30e) \end{align*}$$
		- then with new $a$ and $\psi$
			- $$\Large \begin{align*} & && (30f) \end{align*}$$
			- $$\Large \begin{align*}& \varphi(l,t)\approx\varphi^h(l,t)=\sum_{m=1}^{6N}\psi_m(l)a_m^{\varphi}(t), &&(30f)\\&\delta \varphi(l,t)\approx\delta\varphi^h(l,t)=\sum_{k=1}^{6N}\psi_k(l)b_k^{\varphi}(t), &&(30g) \end{align*}$$
			- for faster computation new formula is
			  background-color:: blue
				- background-color:: blue
				  $$\Large \begin{align*}& \varphi(l,t)\approx\varphi^h(l,t)=\sum_{m=1}^{6}\psi^e_{m}(l)a_{m^*(m,e)}^{\varphi}(t),&&\\&e=\left\{ \begin{align*} &l//\Delta l +1 &&\text{,if}&&l<L \\ &N &&\text{,if}&&l=L \end{align*}\in[1,N] \right. &&\\&m^*(m,e)=m+6(e-1)  &&(30e)  \end{align*}$$
				- for this case:
				  background-color:: blue
					- $$\large \begin{align*} &\psi_i^e(l)=\psi_i^e(\lambda)=\psi_i(\frac{l-\mathcal{N}_{e-1}}{\Delta l}) ,i\in[1,6] && (30b) \end{align*}$$
					- $a_m^{\varphi}$ the same
			-
	- **Boundary Conditions Expressed in the Curve Tangent Angle**
		- Because treating position-based boundary
		  conditions is not directly possible in the curve tangent angle
		  beam model (Eq. 21)
			- this chapter focuses on the
			  development of a strategy to express boundary conditions
			  in terms of higher order derivatives only.
		- Without loss of
		  generality, we consider for our beam model a fixed end,
			- ![image.png](../assets/image_1705604440399_0.png)
		- and, respectively, its dynamic counterpart
			- ![image.png](../assets/image_1705604464721_0.png)
		- Fixing an elastic beam in its position at one end introduces
		  point-wise reaction forces from the mounting onto the beam
		  in x and/or y directions.
		- While these boundary conditions are
		  straight forward to be incorporated in a FEM formulation for
		  a model in the parameters x and y, e.g., Eq. 9a,
			- the FEM
			  description of the reduced model (Eq. 28a) directly acts on
			  the tangent angle function φ and its derivatives
			- thus not
			  offering any parameter to incorporate position boundary
			  conditions.
		- the curve tangent
		  angle dynamics Eq. 21 do not require any position
		  parameters to govern the beam profile.
		- We thus propose to transform the position boundary
		  condition at node N0
			- which cannot be incorporated
			  directly, into a dynamic boundary condition
				- for the neighboring node N 1 entirely expressed in curve tangent
				  coefficients aφ.
		- To define this substitutional boundary
		  condition,
			- we first derive another FEM formulation for
			  the beam model parametrized in Cartesian coordinates
			  (Eq. 9a).
		- Not considering external forces for simplicity
		  and recalling the geometric identities Eq. 11, the beam
		  equation Eq. 9a reduces to
		- ![image.png](../assets/image_1706775871193_0.png)
		- ![image.png](../assets/image_1705605070085_0.png)
		- Formulating the weak forms and applying another integration
		  by parts reads
		- ![image.png](../assets/image_1705605326418_0.png)
		- and the Bubnov–Galerkin approximation for $a_m^{x},\,a_m^{y}$ and this vector exactly the same as $a_m^{\varphi}$ only with $x,\,y$ instead of $\varphi$
			- $$ \begin{align*}& x(l,t)\approx x^h(l,t)=\sum_{m=1}^{6N}\psi_m(l)a_m^{x}(t), &&(35a) && y(l,t)\approx y^h(l,t)=\sum_{m=1}^{6N}\psi_m(l)a_m^{y}(t), &&(35b) \\ &\delta x(l,t)\approx\delta x^h(l,t)=\sum_{k=1}^{6N}\psi_k(l)b_k^{x}(t), &&(35c) &&\delta y(l,t)\approx\delta y^h(l,t)=\sum_{k=1}^{6N}\psi_k(l)b_k^{y}(t), &&(35d) \end{align*}$$
		- with the same set of orthogonal functions ψ that leads to the
		  systems of equations
			- second term is equal zero as
		- $$\Large \begin{align*} &\text{left side x}=\bm{f}^x_1 (\bm{a}^x) + \bm{f}^x_2 (\bm{a}^x) && (36a)\\&\text{left side y}=\bm{f}^y_1 (\bm{a}^y) + \bm{f}^y_2 (\bm{a}^y) && (36b) \end{align*}$$
		- $$\Large \begin{align*} &\text{left side x}=\bm{f}^x_1 (\bm{a}^{\varphi}) + \bm{f}^x_2 (\bm{a}^{\varphi}) && (36.2a)\\&\text{left side y}=\bm{f}^y_1 (\bm{a}^{\varphi}) + \bm{f}^y_2 (\bm{a}^{\varphi}) && (36.2b) \end{align*}$$
		-
		- $$\large \begin{align*} &\bm{f}^x_1 (\bm{a}^x) = -\int_0^L sin(\varphi)(\varphi)_{ll}(\delta x)_l dl && \bm{f}^x_2 (\bm{a}^x)=sin(\varphi)(\varphi)_{ll}\delta x \Bigg|_0^L && (36.1a)\\& \bm{f}^y_1 (\bm{a}^y) =  \int_0^L cos(\varphi)(\varphi)_{ll}(\delta y)_l dl && \bm{f}^y_2 (\bm{a}^y) = -cos(\varphi)(\varphi)_{ll}\delta y \Bigg|_0^L && (36.1b) \end{align*}$$
		-
		- $$\Large \begin{align*} & (a_{m}^{x}(t))_{tt}=\left\{ \begin{align*} &(x)_{tt}(\mathcal{N}_{(m-1)//6},t) &&\text{,if}&&m\%6=1\\           &(x)_{ltt}(\mathcal{N}_{(m-1)//6},t) &&\text{,if}&&m\%6=2\\               &(x)_{lltt}(\mathcal{N}_{(m-1)//6},t) &&\text{,if}&&m\%6=3\\          &(x)_{tt}(\mathcal{N}_{(m-1)//6+1},t) &&\text{,if}&&m\%6=4\\        &(x)_{ltt}(\mathcal{N}_{(m-1)//6+1},t) &&\text{,if}&&m\%6=5\\            &(x)_{lltt}(\mathcal{N}_{(m-1)//6+1},t) &&\text{,if}&&m\%6=0 \end{align*} \right.         ,&&m\in[1,6N] &&  \end{align*}$$
		  $$\Large \begin{align*} &(a_{m}^{y}(t))_{tt}=\left\{ \begin{align*} &(y)_{tt}(\mathcal{N}_{(m-1)//6},t) &&\text{,if}&&m\%6=1\\           &(y)_{ltt}(\mathcal{N}_{(m-1)//6},t) &&\text{,if}&&m\%6=2\\               &(y)_{lltt}(\mathcal{N}_{(m-1)//6},t) &&\text{,if}&&m\%6=3\\          &(y)_{tt}(\mathcal{N}_{(m-1)//6+1},t) &&\text{,if}&&m\%6=4\\        &(y)_{ltt}(\mathcal{N}_{(m-1)//6+1},t) &&\text{,if}&&m\%6=5\\            &(y)_{lltt}(\mathcal{N}_{(m-1)//6+1},t) &&\text{,if}&&m\%6=0 \end{align*} \right.         ,&&m\in[1,6N] && (39a) \end{align*}$$
		- $$\Large \begin{align*} & (a_{m}^{x}(t))_{t}=\left\{ \begin{align*} &(x)_{t}(\mathcal{N}_{(m-1)//6},t) &&\text{,if}&&m\%6=1\\           &(x)_{lt}(\mathcal{N}_{(m-1)//6},t) &&\text{,if}&&m\%6=2\\               &(x)_{llt}(\mathcal{N}_{(m-1)//6},t) &&\text{,if}&&m\%6=3\\          &(x)_{t}(\mathcal{N}_{(m-1)//6+1},t) &&\text{,if}&&m\%6=4\\        &(x)_{lt}(\mathcal{N}_{(m-1)//6+1},t) &&\text{,if}&&m\%6=5\\            &(x)_{llt}(\mathcal{N}_{(m-1)//6+1},t) &&\text{,if}&&m\%6=0 \end{align*} \right.         ,&&m\in[1,6N] &&  \end{align*}$$
		  $$\Large \begin{align*} &(a_{m}^{y}(t))_{t}=\left\{ \begin{align*} &(y)_{t}(\mathcal{N}_{(m-1)//6},t) &&\text{,if}&&m\%6=1\\           &(y)_{lt}(\mathcal{N}_{(m-1)//6},t) &&\text{,if}&&m\%6=2\\               &(y)_{llt}(\mathcal{N}_{(m-1)//6},t) &&\text{,if}&&m\%6=3\\          &(y)_{t}(\mathcal{N}_{(m-1)//6+1},t) &&\text{,if}&&m\%6=4\\        &(y)_{lt}(\mathcal{N}_{(m-1)//6+1},t) &&\text{,if}&&m\%6=5\\            &(y)_{llt}(\mathcal{N}_{(m-1)//6+1},t) &&\text{,if}&&m\%6=0 \end{align*} \right.         ,&&m\in[1,6N] && (39b) \end{align*}$$
		- While the right hand side of the FEM formulations (Eq. 36) is already defined in the curve tangent angle φ, what remains is to also rewrite coefficient vectors in terms of φ instead of x and y. Starting again from the geometric identities (x)l ≡ cos(φ) and (y)l ≡ sin(φ), the coefficients in Eq. 39 can be expressed as
		- $$ \begin{align*} &x(l,t)=x(0,t)+\int_0^l cos(\varphi(s,t))ds && (40.1a)\\&y(l,t)=y(0,t)+\int_0^l sin(\varphi(s,t))ds && (40.1b) \end{align*}$$
		- $$ \begin{align*} &(x)_{t}(l,t)=(x)_{t}(0,t)-\int_0^l sin(\varphi(s,t))(\varphi(s,t))_{t}ds && (40.1.1a)\\&(y)_{t}(l,t)=(y)_{t}(0,t)+\int_0^l cos(\varphi(s,t))(\varphi(s,t))_{t}ds && (40.1.1b) \end{align*}$$
		- $$ \begin{align*} &(x)_{tt}(l,t)=(x)_{tt}(0,t)-\int_0^l cos(\varphi(s,t))(\varphi(s,t))^2_{t}ds-\int_0^l sin(\varphi(s,t))(\varphi(s,t))_{tt}ds && (40.1.2a)\\&(y)_{tt}(l,t)=(y)_{tt}(0,t)-\int_0^l sin(\varphi(s,t))(\varphi(s,t))^2_{t}ds+\int_0^l cos(\varphi(s,t))(\varphi(s,t))_{tt}ds && (40.1.2b) \end{align*}$$
		- $$ \begin{align*} &(x(l,t))_{l}= cos(\varphi(l,t)) && (40.2a)\\&(y(l,t))_{l}=sin(\varphi(l,t)) && (40.2b) \end{align*}$$
		- $$ \begin{align*} &(x(l,t))_{lt}= -sin(\varphi(l,t))(\varphi(l,t))_{t} && (40.2.1a)\\&(y(l,t))_{lt}=cos(\varphi(l,t))(\varphi(l,t))_{t} && (40.2.1b) \end{align*}$$
		- $$ \begin{align*} &(x(l,t))_{ltt}= -cos(\varphi(l,t))(\varphi(l,t))^2_{t}-sin(\varphi(l,t))(\varphi(l,t))_{tt} && (40.2.2a)\\&(y(l,t))_{ltt}=-sin(\varphi(l,t))(\varphi(l,t))^2_{t}+cos(\varphi(l,t))(\varphi(l,t))_{tt} && (40.2.2b) \end{align*}$$
		- $$ \begin{align*} &(x(l,t))_{ll}= -sin(\varphi(l,t))(\varphi(l,t))_{l} && (40.3a)\\&(y(l,t))_{ll}=cos(\varphi(l,t))(\varphi(l,t))_{l} && (40.3b) \end{align*}$$
		- $$ \begin{align*} &(x(l,t))_{llt}= -cos(\varphi(l,t))(\varphi(l,t))_{t}(\varphi(l,t))_{l}-sin(\varphi(l,t))(\varphi(l,t))_{lt} && (40.3.1a)\\&(y(l,t))_{llt}=-sin(\varphi(l,t))(\varphi(l,t))_{t}(\varphi(l,t))_{l}+cos(\varphi(l,t))(\varphi(l,t))_{lt} && (40.3.1b) \end{align*}$$
		- $$ \begin{align*} &(x(l,t))_{lltt}=  sin(\varphi(l,t))(\varphi(l,t))^2_{t}(\varphi(l,t))_{l}-&&\\&-cos(\varphi(l,t))(\varphi(l,t))_{tt}(\varphi(l,t))_{l} - cos(\varphi(l,t))(\varphi(l,t))_{t}(\varphi(l,t))_{lt}- &&\\& -cos(\varphi(l,t))(\varphi(l,t))_{t}(\varphi(l,t))_{lt}-sin(\varphi(l,t))(\varphi(l,t))_{ltt} && (40.3.2a)\\&(y(l,t))_{lltt}=-cos(\varphi(l,t))(\varphi(l,t))^2_{t}(\varphi(l,t))_{l}-&&\\&-sin(\varphi(l,t))(\varphi(l,t))_{tt}(\varphi(l,t))_{l}-sin(\varphi(l,t))(\varphi(l,t))_{t}(\varphi(l,t))_{lt}+&&\\&-sin(\varphi(l,t))(\varphi(l,t))_{lt}(\varphi(l,t))_{t}+cos(\varphi(l,t))(\varphi(l,t))_{ltt} && (40.3.2b) \end{align*}$$
		- 9 equations got containing Volterra integrals with an upper limit l.
		- Recalling the spline approximation φh from Eq. 30fg and using it for the acceleration terms,
		- Eventually, Eq. 38 can be evaluated entirely in φ with the left hand sides
		- ![image.png](../assets/image_1705681246173_0.png)
		- ![image.png](../assets/image_1705681265531_0.png)
		- ![image.png](../assets/image_1705681227558_0.png)
-
-
