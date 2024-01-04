package progistar.pXg.data;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Hashtable;

import org.ahocorasick.trie.Emit;
import org.ahocorasick.trie.Trie;

public class Fasta {

	private ArrayList<FastaRecord> records = new ArrayList<FastaRecord>();
	
	public Fasta (String fileName) {
		if(fileName == null) {
			return;
		}
		
		File file = new File(fileName);
		
		// check if the file is available. 
		if(file == null || !file.exists()) {
			return;
		}
		
		try {
			boolean isWarning = false;
			BufferedReader BR = new BufferedReader(new FileReader(file));
			String line = null;
			
			FastaRecord fastaRecord = null;
			StringBuilder sequence = new StringBuilder();
			while((line = BR.readLine()) != null) {
				if(line.startsWith(">")) {
					// header parsing
					if(fastaRecord != null) {
						fastaRecord.sequence = sequence.toString();
					}
					fastaRecord = new FastaRecord();
					String id = line.split("\\s")[0].substring(1);;
					
					// check if the id contains "|"
					if(id.contains("|")) {
						if(!isWarning) {
							System.out.println("Detected | character in a header of the fasta file");
							System.out.println("The character is converted :");
							isWarning = true;
						}
						id = id.replace("|", ":");
					}
					
					fastaRecord.id = id;
					records.add(fastaRecord);
					sequence.setLength(0);
				} else {
					sequence.append(line.replace("I", "L"));
				}
			}
			
			if(fastaRecord != null && fastaRecord.sequence == null) {
				fastaRecord.sequence = sequence.toString();
			}
			
			BR.close();
		}catch(IOException ioe) {
			
		}
	}
	
	/**
	 * Find all fasta records matching to the given peptide sequences. <br>
	 * Return Hashtable consisting with:<br>
	 * key: peptide seuqence <br>
	 * value: ArrayList of matched IDs <br>
	 * 
	 * 
	 * @param sequences
	 * @return
	 */
	public Hashtable<String, ArrayList<String>> findAll (ArrayList<String> sequences) {
		Hashtable<String, ArrayList<String>> matchedList = new Hashtable<String, ArrayList<String>>();
		
		Trie trie = Trie.builder().addKeywords(sequences).build();
		
		// find all fasta records
		this.records.forEach(record -> {
			Collection<Emit> emits = trie.parseText(record.sequence);
			
			for(Emit emit : emits) {
				String pSeq = emit.getKeyword();
				
				ArrayList<String> ids = matchedList.get(pSeq);
				if(ids == null) {
					ids = new ArrayList<String>();
					matchedList.put(pSeq, ids);
				}
				ids.add(record.id);
			}
		});
		
		return matchedList;
	}
}
