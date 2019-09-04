filename = 
"/localdisk/home/PATH_TO_FOLDER/fasta_file_for_subtyping.fasta";

inArguments 
= {"0":"Nucleotide",  
   "1":filename, 
   "2":filename+".out", 
   "3":"Summary+Detail",
   "4":filename+ ".details"};

ExecuteAFile ("/localdisk/software/SCUEAL_fullG/SCUEAL-master/TopLevel/MPIScreenFASTA.bf",inArguments);
