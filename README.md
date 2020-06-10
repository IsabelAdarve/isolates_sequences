# isolates_sequences
Here you will find the Python script used for generating the trimmed sequences for every .ab1 file of the DNA samples that were sequenced by Sanger method in the Adarve-Rengifo et al., 2020.

After obtaining the results of the sequencing by the Sanger method, folders with four files for every sequenced DNA sample of extentions .ab1 , .pdf , .phd and .txt are sent to the investigator by the laboratory in charge of doing the sequencing. 

To proceed with the identification of every sample:

First a work folder must be created containing the folder with data sent by the laboratory and the python script. Then from the terminal using the cd command move to the work folder where the information generated by the analysis script will be saved. The following is the syntax to execute the python script to process Sanger data in order to obtain trimmed sequences 
$ python readmacrogene.py outputfolder 6 30. Every parameter stands for:

  * readmacrogene.py: the python script will be using.
  * outputfolder: folder sent by the laboratory with data sent after sequencing. 
  * 6: window of 6 bp (this value can be modified by who's running the code).
  * 30: score for Phred value (this value can be modified by who's running the code).
  
Files are listed to extract only the .ab1 files from where the DNA sequences will be obtained for later timming. Posterior to files extraction and their respective trim, these new generated files will be .txt and saved into 'results'* independent folder. Final files will be able to be uploaded into the Basic Local Alignment Search Tool (BLAST) software at National Center for Biotechnology Information (NCBI) for the taxonomic characterization and the organism with complete genome and higher 'Total Score' will be selected (these parameters are freely selectable according to the user).

*NOTE: for this version the folder 'results' must previously exist on the work folder.
