Directory containing all workflow for bridge membership study.

~/superimposition/superimposition.py
------------------------------------

This script iterates over all files in the PDB and attempts
to compare the relationship between 2-bridges of form
A : Met : B and closely spaced redox active aromatic chains.
A, B = any of Phe, Tyr Trp.

(1) Import PDB file from PDB
(2) Search for bridges and search for chains
(3) Isolate 2-bridges
(4) Does PDB file contain both 2-bridge and chain?
    if YES;
      continue
    else
      next PDB file
(5) Assess relationship between chain(s) and bridge(s)
(6) Output relationship to text file

Example entry in text file looks as follows:
===============================================================
1A8D : Chains  : TYR374, TYR348,  
1A8D : Chains  : TYR287, TYR285,  
1A8D : Chains  : TRP168, TYR388, TYR245, TYR242, TRP92,  
1A8D : Chains  : TYR266, TRP449, TRP255, TYR307, TYR303, TYR294,  
1A8D : Chains  : TYR427, TRP426,  
1A8D : Chains  : TRP139, TRP128,  
1A8D : Bridges : TYR80,  PHE84,  
1A8D : Bridges : TYR163, PHE153,  
===============================================================


~/superimposition/analyze.py
----------------------------

This script takes all files of form:
===============================================================
1A8D : Chains  : TYR374, TYR348,  
1A8D : Chains  : TYR287, TYR285,  
1A8D : Chains  : TRP168, TYR388, TYR245, TYR242, TRP92,  
1A8D : Chains  : TYR266, TRP449, TRP255, TYR307, TYR303, TYR294,  
1A8D : Chains  : TYR427, TRP426,  
1A8D : Chains  : TRP139, TRP128,  
1A8D : Bridges : TYR80,  PHE84,  
1A8D : Bridges : TYR163, PHE153,  
===============================================================

And compares the relationship between all bridges and all chains:

"NR : {' PHE84', ' TYR80'} : {' TYR348', ' TYR374'}",
"NR : {' PHE84', ' TYR80'} : {' TYR285', ' TYR287'}",
"NR : {' PHE84', ' TYR80'} : {' TRP168', ' TRP92', ' TYR388', ' TYR242', ' TYR245'}",
"NR : {' PHE84', ' TYR80'} : {' TRP449', ' TYR266', ' TYR303', ' TRP255', ' TYR307', ' TYR294'}",
"NR : {' PHE84', ' TYR80'} : {' TYR427', ' TRP426'}",
"NR : {' PHE84', ' TYR80'} : {' TRP139', ' TRP128'}",
"NR : {' TYR163', ' PHE153'} : {' TYR348', ' TYR374'}",
"NR : {' TYR163', ' PHE153'} : {' TYR285', ' TYR287'}",
"NR : {' TYR163', ' PHE153'} : {' TRP168', ' TRP92', ' TYR388', ' TYR242', ' TYR245'}",
"NR : {' TYR163', ' PHE153'} : {' TRP449', ' TYR266', ' TYR303', ' TRP255', ' TYR307', ' TYR294'}",
"NR : {' TYR163', ' PHE153'} : {' TYR427', ' TRP426'}",
"NR : {' TYR163', ' PHE153'} : {' TRP139', ' TRP128'}"
 ^
 |
 Then count along this column and export data to pie chart.
 
