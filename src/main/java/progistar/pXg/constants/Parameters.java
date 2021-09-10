package progistar.pXg.constants;

public class Parameters {
	// Input file paths
	public static final String gAnnPath				=	"gann".toLowerCase();
	public static String genomicAnnotationFilePath	=	"C:\\Bioinformatics\\0.Databases\\1.HumanRNAs\\gencode.v19.annotation.sorted.gtf";
	
	public static final String gSeqPath				=	"gseq".toLowerCase();
	public static String sequenceFilePath			=	"C:\\Users\\progi\\Desktop\\Projects\\pXg\\chr1toy.sam";
	
	public static final String pAnnPath				=	"pann".toLowerCase();
	public static String peptideFilePath			=	"";
	
	// Output file path
	public static final String oPath				=	"output".toLowerCase();
	public static String outputFilePaht				=	"";
	
	// GTF partition size
	public static final String gPartitionSize		=	"gann_size".toLowerCase();
	public static int partitionSize					=	5000000; // 5 * 10^6 * 10 * 8 = 0.4 G
	
	public static final String gSeqReadSize			=	"gseq_size".toLowerCase();
	public static int readSize						=	1000000; // 1 * 10^7
	
	// READ sequencing 
	public static byte READ_STRAND					=	Constants.FORWARD_STRAND_READS;
	
	
	// Peptide file
	// for user-friendly purpose, peptideColumnIndex is taken one-based and converted to zero-based.
	public static final String pColumnIndex			=	"pcol".toLowerCase();
	public static int peptideColumnIndex			=	0; // user-specific peptide index
	
	public static final String pParserRegExr		=	"aareg".toLowerCase();
	public static String peptideParserRegExr		=	"[A-Z]"; // read sequence matched to the RegExr.
	// Note that
	// comment is different from field.
	// In pXg definition, field indicates column names and comment shows some meta-data.
	// Comment line must be present on the top of records.
	// Field must follow "comment" when the comment presents in the file.
	public static final String cMarker				=	"cm".toLowerCase();
	public static String commentMarker				=	"#|@|%"; // if line starts with the pattern, the line will be skipped during parsing the file.
	
	// System Parameters
	public static final String numOfThreads			=	"num_threads".toLowerCase();
	public static int nThreads						=	1;
	
	// BLAST PARAMETER
	public static final String makeblastdbPath		=	"makedb".toLowerCase();
	public static String TBLASTN_EXE_PATH			=	"C:\\Bioinformatics\\0.utils\\blast-2.9.0+\\bin\\tblastn.exe";
	
	public static final String tblastnPath			=	"tblastn".toLowerCase();
	public static String MAKEDB_EXE_PATH 			=	"C:\\Bioinformatics\\0.utils\\blast-2.9.0+\\bin\\makeblastdb.exe";
	
}
