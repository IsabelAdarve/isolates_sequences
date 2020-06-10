# isolates_sequences
Here you will find the Python script used for generating the trimmed sequences for every .ab1 file of the DNA samples that were sequenced by Sanger method in the Adarve-Rengifo et al., 2020.

After obtaining the results of the sequencing by the Sanger method, folders with four files for every sequenced DNA sample of extentions .ab1 , .pdf , .phd and .txt are sent

To proceed with the identification of every sample:

From the terminal using the cd command move to the folder were you desire to save the information. Call python3, the script, folder with data sent after sequencing, the desired 'n' and the 'score': folders are listed to extract only the .ab1 files from where the DNA sequences will be obtained for later timming using a window of 6 pdb and a 30 score (these values can be modified by who's running the code).

Posterior to files extraction and their respective trim, these new generated files will be .txt and saved into 'results'* independent folder. Final files will be able to be uploaded into the Basic Local Alignment Search Tool (BLAST) software at National Center for Biotechnology Information (NCBI) for the taxonomic characterization and the organism with complete genome and higher 'Total Score' will be selected (these parameters are freely selectable according to the user).

* NOTE: for this version the folder 'results' must previously exist on the folder were data folder and the script are.
