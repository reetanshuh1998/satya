# tagged-neutron-DIS
Event generator for leading neutron tagged DIS process, which is useful for constraining pion structure function

# 1. Event generation of e p --> e' n X data

**Using ROOT to execute the codes:**
>root -l -q test.cpp

**Using ACLiC to compile and execute:**
>root -l

[0] .x test.cpp+

**Using g++ to compile and execute:**
>sed 's/test()/main()/' test.cpp > test2.cpp

>g++ test2.cpp -o run_sim \`root-config --cflags\` \`root-config --libs\`

>./run_sim


# 2. Get the pionic parton distributions from IMParton

<testPiPDF.cpp> is the sample code for the users.
Using ROOT to execute the codes:                                                                                                                                      
>root -l testPiPDF.cpp


# Reference

Please check G. XIE et al., "Simulation of neutron-tagged deep inelastic scattering at EicC",
Chinese Physics C 45 (2021) 5, 053002 [arXiv:2009.04956],
https://inspirehep.net/literature/1816055
for more information about the event generator.


Download ROOT and get started: https://root.cern.ch/

coding issue contact: rwang@impcas.ac.cn

