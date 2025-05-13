# ProtAnnoScripts
This repository contains Python functions to parse the output from various protein annotation tools and transform them into pandas dataframes.  
They are used for protein product prediction, protein family classification, functional annotation submission to GenBank and prediction of subcellular localization.   



## Subcellular localization
The parse_localization.annotate_location function predicts the localization of proteins by integrating annotations from UniFire, DeepLoc2, SignalP6, TMHMM2 and NetGPI. Settings are optimized for fungal proteins the following way:

1. SignalP6 && NetGPI ->  Extracellular GPI-anchored
2. DeepLoc2['Extracellular'] > 0.52753906 -> Extracellular
3. UniFire -> Intracellular, Extracellular, Membrane (optional)
4. SignalP6 && noTM -> Extracellular
5. TMHMM | DeepLoc2['Membrane'] > 0.64638672 -> Membrane
6. Deeploc(Nucles|Cytoplasm) -> Intracellular
7. Rest remains blank but considered Intracellular


## Pending to explain
- annie-style file generation
- classification in enzyme groups
- CAZyme annotations