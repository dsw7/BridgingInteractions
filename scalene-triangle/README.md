### All files for the scalene triangle study  
---  
### Description of directories:  

    ~/scalene-triangle/pymol-get-surface-example // example PyMOL API code depicting how surface coordinates were obtained
    ~/scalene-triangle/coordinate-to-surface     // example code depicting how surface to MT or SD face was chosen
    ~/scalene-triangle/libs                      // contains all dependencies
    ~/scalene-triangle/legacy                    // some old and "hacky" code from the Summer 2018 superimposition study
    ~/scalene-triangle/scalene_main.py           // the "main" script that performs the superimposition analysis 
                                                 // for PDB entries located in delim.txt
    ~/scalene-triangle/delim.txt                 // a .txt file that tells main script which entries to analyze
    ~/scalene-triangle/analyze.py                // the main script loads data into a .txt file. This script analyzes
                                                 // that .txt file and generates three sets of histograms
    ~/scalene-triangle/analyze_sidebyside.py     // almost identical to analyze.py but subplots are rendered side-by-side
                                                 // plot dimensions, font sizes & tick frequency also changed                                                                                           
    ~/scalene-triangle/test.py                   // compares data from Summer 2018 study to Winter 2018 study    
    ~/scalene-triangle/results_scalene.txt       // the source file for analyze.py / analyze_sidebyside.py scripts             
    
---  
### analyze_sidebyside.py output examples:  
  
<img src="https://github.com/dsw7/BridgingInteractions/blob/master/scalene-triangle/all_phe.png">  
<img src="https://github.com/dsw7/BridgingInteractions/blob/master/scalene-triangle/phe_aro.png">  
<img src="https://github.com/dsw7/BridgingInteractions/blob/master/scalene-triangle/aro_aro.png">  
