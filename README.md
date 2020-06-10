# isolates_sequences
Here you will find the Python script used for generating the trimmed sequences for every .ab1 file of the DNA samples that were sequenced by Sanger method in the Adarve-Rengifo et al., 2020.

After obtaining the results of the sequencing by the Sanger method three folders were sent (each one with four files for every sequenced isolate) of extentions .ab1 , .pdf , .phd and .txt

To proceed with the identification of every isolate:

From the terminal using the cd command move to folder 'finalCompu2'. Call python3, the script, 'data' folder, the desired 'n' and the 'score': all 3 folders were listed to extract only the .ab1 files from where the DNA sequences will be obtained for later timming using a window of 6 pdb and a 30 score (these values can be modified for who's running the code).

Posterior to files extraction and their respective trim, these new generated files will be .txt and saved into 'results'* independent folder. Final files will be able to be uploaded into the Basic Local Alignment Search Tool (BLAST) software at National Center for Biotechnology Information (NCBI) for the taxonomic characterization and the organism with complete genome and higher 'Total Score' will be selected (these parameters are freely selectable according to the user).

* NOTE: for this version the folder 'results' must previously exist on the folder were data folder and the script are.
