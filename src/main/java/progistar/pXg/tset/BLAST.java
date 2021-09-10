package progistar.pXg.tset;

import java.io.IOException;

import progistar.pXg.constants.Parameters;

public class BLAST {

	// OUTFMT 6 format INDEX
	public static final int QUERY_INDEX = 0;
	public static final int SUBJECT_INDEX = 1;
	public static final int PIDENT_INDEX = 2;
	public static final int ALIGNMENT_LENGTH_INDEX = 3;
	public static final int MISMATCH_INDEX = 4;
	public static final int GAP_INDEX = 5;
	public static final int QUERY_START_INDEX = 6;
	public static final int QUERY_END_INDEX = 7;
	public static final int SUBJECT_START_INDEX = 8;
	public static final int SUBJECT_END_INDEX = 9;
	public static final int EVALUE_INDEX = 10;
	public static final int BIT_SCORE_INDEX = 11;
	
	/**
	 * inputFasta path should be absolute path. <br>
	 * dbtype: nucl | prot <br>
	 * 
	 * @param inputFasta
	 * @param dbtype
	 * @throws IOException
	 */
	public static void makeBlastDB (String inputFasta, String dbtype) {
		String makeBlastExe = Parameters.MAKEDB_EXE_PATH;
		
		try {
			String cmd = makeBlastExe+" -in "+ inputFasta +" -input_type fasta -dbtype "+dbtype;
			Runtime runtime = Runtime.getRuntime();
			Process makeDB = runtime.exec(cmd);
			System.out.println(cmd);
			makeDB.waitFor();
		}catch(Exception e) {}
	}
	
	/**
	 * The query format must be Fasta. <br>
	 * 
	 * @param fasta
	 * @param dbPath
	 */
	public static void tblastn (String fastaQuery, String dbPath, String outputPath, boolean isOutfmt6) {
		String outfmt6 = isOutfmt6 ? "-outfmt 6" : "";
		
		// write query into tempFile
		try {
			// blastn
			String tblastnExe = Parameters.TBLASTN_EXE_PATH;
			// short sequence:
			// evalue is really important.!
			// PAM30 is the most used matrix for short-sequence
			String cmd = tblastnExe+" -query "+ fastaQuery +" -db "+dbPath +" -evalue 1000 -word_size 3 -num_alignments 1000 -max_hsps 100 "+outfmt6+" > " + outputPath ;
			Runtime runtime = Runtime.getRuntime();
			Process makeDB = runtime.exec(cmd);
			System.out.println(cmd);
			makeDB.waitFor();
		}catch(Exception e) {}
	}
}
