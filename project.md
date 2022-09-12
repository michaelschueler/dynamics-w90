src_dir: ./src
output_dir: ./doc
exclude: Mlebedev_weights.F90
exclude_dir: ./src/math/lebedev/
author: Michael Schueler
author_description: The `dynamics-w90` package: Berry phase properties, light-induced dynamics, photoemission & more from Wannier functions.

Hello ${USER}, welcome to the `dynamics-w90` package! The package contains a number of programs to compute various observables from Wannier functions, which are read output from (Wannier90)[http://www.wannier.org]. 

At the moment `dynamics-w90` contains the following main programs:

1. [[wann_evol]] / [[wann_evol_mpi]]: Computes the time-dependent dynamics of an electron system described by a Wannier Hamiltonian upon laser excitation and computes various observables. Direct implementation of the gauge-invariant formulation from [Phys. Rev. B 103, 1155409 (2021)](https://link.aps.org/doi/10.1103/PhysRevB.103.155409) 
2. [[wann_calc]]: Calculation of band properties including orbital weights, spin texture, Berry curvature and orbital angular momentum (OAM). We implemented the modern theory of OAM to include all non-local contributions.
3. [[wann_prune]]: Compresses a Wannier Hamiltonian by cutting small hopping amplitudes. Output can be written in original Wannier90 or in hdf5 format.
4. [[wann_soc]]: Extends a Wannier Hamiltonian in spin space and adds atomic spin-orbit coupling (SOC).
5. [[arpes]] / [[arpes_mpi]]: Computes angle-resolved photoemission spectrum from Wannier functions. 


@Note The code is under active development. Other functionalities will be added soon!

@Bug You can have multi-paragraph versions of these too! That means you can include

ordered lists
unordered lists
images
etc.
Isn't that cool? @endbug

@Bug Hey I'm doing it again...

This ones ends mid...@endbug ...paragraph.

You can have as many paragraphs as you like here and can use
headlines, links, images, etc. Basically, you can use anything in
Markdown and Markdown-Extra. Furthermore, you can insert LaTeX into
your documentation. So, for example, you can provide inline math using
like ( y = x^2 ) or math on its own line like \(x = \sqrt{y}\) or $$ e
= mc^2. $$ You can even use LaTeX environments! So you can get
numbered equations like this:
\begin{equation}
\left[i \partial_t - h(t) \right] G(t,t^\prime)  = \int_{C} d \bar{t} \Sigma(t,\bar{t})G(\bar{t},t^\prime)
\end{equation}
So let your imagination run wild. As you can tell, I'm more or less
just filling in space now. This will be the last sentence.
