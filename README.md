#Molecular Draw


<img width="1920" height="1080" alt="Screenshot 2026-01-01 at 09 35 22 (2)" src="https://github.com/user-attachments/assets/e84f183e-1f69-40ba-86c2-cf552368cc0d" />


I built a system that converts a textual chemical description into a physically realistic, interactive 3D molecular model that runs end-to-end from backend computation to frontend visualization. The project started from studying how chemists represent molecules using **SMILES notation**, and how that abstract string encodes real chemical structure, bonding, and reactivity. While working on this, I also studied **organic chemistry concepts** such as **addition reactions, esterification, functional groups, bond formation, and steric effects**, and mapped those ideas directly into computational representations.

On the backend, I implemented a molecular processing service using **Python and RDKit**. When a SMILES string such as `"CCO"` is submitted, the system parses it into a molecular graph where atoms are nodes and bonds are edges. This mirrors how chemical connectivity is reasoned about during reactions, where the arrangement of atoms and bond orders determines how molecules transform. Since SMILES omits hydrogens by default, I explicitly added hydrogen atoms to satisfy valency rules, which is essential for chemically correct geometry and force-field calculations. This step directly reinforced concepts like carbon tetravalency and hydrogen participation in reactions such as esterification and nucleophilic addition.

To move from a purely topological structure to a real spatial molecule, I generated **3D atomic coordinates** using RDKit’s distance geometry embedding. The algorithm constructs a distance matrix based on known bond lengths, angles, and stereochemical constraints, then places atoms in three-dimensional space while satisfying those constraints. This allowed me to study how **molecular shape, bond angles, and spatial orientation** influence reactivity, particularly steric hindrance and functional group accessibility.

Once the initial geometry was generated, I applied **MMFF94 force-field optimization** to minimize the molecule’s energy. This step reduces unrealistic bond strain, torsional strain, and van der Waals clashes, producing a stable conformation that reflects how real molecules settle into low-energy states. Working with force-field optimization helped me connect chemical theory—such as why molecules adopt certain shapes—to computational energy models.

After optimization, I extracted the final atomic coordinates and bond information and converted them into a structured JSON format. Each atom includes its element type and precise x, y, z coordinates (in angstroms), and each bond records its connectivity and bond order. This clean data model allowed the backend to remain chemistry-focused while exposing a simple API for visualization. I deployed the backend as a service, ensuring deterministic outputs, reproducibility, and efficient response times.

On the frontend, I built an interactive molecular viewer using **React, Three.js, and React Three Fiber**. The 3D scene is rendered using WebGL with physically based lighting and a perspective camera. Atoms are visualized as spheres positioned directly from the backend-computed coordinates, colored using the **CPK color scheme** so chemical identity is immediately recognizable. Atomic radii are derived from element-specific or van der Waals radii, making spatial size and crowding visually apparent.

Bonds are rendered as cylinders connecting atoms, with orientation calculated using vector math and quaternions to align each bond correctly in 3D space. Bond thickness and appearance reflect bond order, allowing clear visual distinction between single, double, and triple bonds. This made it easy to visually inspect bonding patterns and reason about reaction pathways and bond formation.

I implemented multiple molecular representation modes—ball-and-stick, wireframe, stick, van der Waals, and line models—each emphasizing different chemical properties. For example, van der Waals representations highlight steric bulk and molecular packing, while ball-and-stick models emphasize connectivity and functional groups. These modes directly supported my study of how molecular structure affects chemical behavior.

To make the visualization usable and exploratory, I added interactive orbit controls for rotation, zooming, and panning, along with atom selection and highlighting. Clicking an atom updates application state, enabling future extensions such as structure editing, reaction visualization, or functional group substitution.

Overall, this project allowed me to deeply study and implement:

* Chemical structure representation using SMILES
* Atomic bonding, valency, and hydrogen saturation
* 3D molecular geometry and conformational stability
* Force-field energy minimization (MMFF94)
* Functional groups and reaction concepts like **addition reactions and esterification**
* The relationship between chemical theory and computational modeling
* Backend API design and deployment for scientific computation
* Real-time 3D rendering and interaction in the browser

The final result is a deployed, full-stack system that transforms a simple text string into a chemically accurate, interactive 3D molecule, allowing users to visually explore structure, bonding, and spatial relationships in real time.
