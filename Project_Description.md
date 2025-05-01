# MasterThesis
It contains the scripts made during my Thesis at LACDR, working on computational models of phage-antibiotic combination therapy in case of collateral effects.
Collateral effects are evolutionary phenomena where gaining resistance to one treatment changes the sensitivity to the other.
In this case we will focus on collateral sensitivity (CS) caused by phage resistance, i.e. increased sensitivity to the antibiotic while bacteria become more sensitive to phages.
The model was developed in R using the DeSolve package and the RxOde2 package. We are in a in vitro scenario where we have 4 bacterial strains, phages and the antibiotic. The 4 strains are the "wild type" strain, antibiotic resistant strain, phage resistant strain and double resistant strain.
We first evaluate simultaneous administration, we then focus on sequential administration, varying the time delay between treatments and order of treatment. We then propose a 3 day treatment scheme based on previous findings and finally evaluate the fitness costs of resistance and mutation rates values in the setting of simultaneous administration.
