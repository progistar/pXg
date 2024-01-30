package progistar.pXg.constants;

public class Constants {
	// [major].[minor].[patch].
	// major: significant change. it strongly affects to performance.
	//        plus, basically, this version number intends to upgrade behaviors of tool.
	// minor: critical bug fix, change in results something like that.
	// patch: tiny bug fix, typo, and change code styles something like that.
	public static final String VERSION = "pXg v2.2.0";
	public static final String RELEASE = "(release 2024-01-27)";
	public static final String INTRODUCE = "Seunghyuk Choi and Eunok Paek in the Department of Computer Science at Hanyang University in Seoul, South Korea. ";
	
	public static long	UNIQUE_RUN_ID	=	0;
	
	// Translation method
	public static final int THREE_FRAME	=	3;
	public static final int SIX_FRAME	=	6;
	
	// Decoy prefix
	public static final String DECOY_PREFIX	=	"Decoy$";

	//
	public static final byte EXON = 30;
	
	// EXON is classified into three types such as CDS, UTR and NCDS (Non CoDing Sequence)
	public static final byte CDS = 100;
	public static final byte UTR5 = 25;
	public static final byte UTR3 = 23;
	public static final byte NCDS = 20;
	public static final byte INTRON = 0;
	public static final byte INTERGENIC = -4;
	
	// Regional Character
	public static final char MARK_CDS 			=	'C';
	public static final char MARK_5UTR			=	'F';
	public static final char MARK_3UTR			=	'T';
	public static final char MARK_NCDS			=	'N';
	public static final char MARK_INTRON		=	'I';
	public static final char MARK_INTERGENIC	=	'-';
	
	public static final char MARK_MAPPED		=	'M';
	public static final char MARK_SOFTCLIP		=	'?';
	public static final char MARK_UNMAPPED		=	'*';
	
	// Alternative Splicing Character
	public static final char MARK_AS			=	'A'; // for alternative splicing form
	public static final char MARK_CA			=	'C'; // for canonical form
	
	// Transcript coding type
	public static final byte NON_CODING_TRANSCRIPT	=	0;
	public static final byte CODING_TRANSCRIPT		=	1;
	
	// Frame type
	public static final char NO_FRAME		=	'N';
	public static final char IN_FRAME		=	'I';
	public static final char OUT_OF_FRAME	=	'O';
	
	// Frame index
	// Note that FRAME_X denots NO_FRAME.
	public static final byte FRAME_0		=	0;
	public static final byte FARME_1		=	1;
	public static final byte FRAME_2		=	2;
	public static final byte FRAME_X		=	3;
	
	// Mutation type
	public static final byte SNP			=	0;
	public static final byte INS			=	1;
	public static final byte DEL			=	2;
	
	// Mutation status
	public static final String MUTATION_MISSENSE		=	"missense";
	public static final String MUTATION_SILENT	=	"silent";
	public static final String MUTATION_STOPLOSS		=	"stoploss";
	public static final String MUTATION_INDELS			=	"indels";
	public static final String MUTATION_STOPGAIN		=	"nonsense"; // currently, it is not supported.
	
	
	// RNA-Seq parameter
	public static final byte FORWARD_STRAND_READS	=	0;
	public static final byte REVERSE_STRAND_READS	=	1;
	
	// RNA-Seq Mock
	public static final byte MOCK_NONE			=	0;
	public static final byte MOCK_REVERSE		=	1;
	public static final byte MOCK_PSD_REVERSE		=	2;
	
	public static final byte MOCK_ALL			=	0; // use all multiply-assigned read count
	public static final byte MOCK_MAX_ONE		=	1; // use maximum one as a representative
	public static final byte MOCK_MEAN			=	2; // use mean as a representative
	
	// TASKS
	public static final int TASK_G_MAP				=	1;
	
	// Output Annotation
	public static final String OUTPUT_G_UNIQUE_ID	=	"[ID]";
	public static final String OUTPUT_G_SEQUENCE	=	"[SEQ]";
	public static final String OUTPUT_G_REGION		=	"[REGION]";
	public static final String OUTPUT_G_PEPTIDE		=	"[PEPTIDE]";
	public static final String OUTPUT_G_QSCORE		=	"[QSCORE]";
	
	// Events
	public static final String EVENT_ANTISENSE		=	"asRNA";
	public static final String EVENT_5UTR			=	"5`-UTR";
	public static final String EVENT_3UTR			=	"3`-UTR";
	public static final String EVENT_NONCODING		=	"ncRNA";
	public static final String EVENT_INTERGENIC		=	"IGR";
	public static final String EVENT_INTRON			=	"IR";
	public static final String EVENT_FRAMESHIFT		=	"FS";
	public static final String EVENT_PROTEINCODING	=	"PC";
	public static final String EVENT_AS				=	"AS";
	public static final String EVENT_UNKNOWN		=	"unknown";
	public static final String EVENT_SOFTCLIP		=	"Softclip";
	
	// PSM Status
	public static final byte PSM_STATUS_UNDEF		=	-1;
	public static final byte PSM_STATUS_DECOY		=	0;
	public static final byte PSM_STATUS_BOTH		=	1;
	public static final byte PSM_STATUS_TARGET		=	2;
	
	// ID status
	public static final String ID_NULL				=	"-";
	
	
	/** pXg TSV Parser Constants **/
	public static final String INFERRED_PEPTIDE_COLUMN_NAME	=	"InferredPeptide";
	public static final String NUCLEOTIDE_COLUMN_NAME		=	"ObservedNucleotide";
	public static final String STRAND_COLUMN_NAME			=	"Strand";
	public static final String EVENT_COLUMN_NAME			=	"Events";
	
	// hidden parameters for revision
	public static final String CAL_PHRED_MIN				=	"min";
	public static final String CAL_PHRED_MAX				=	"max";
	public static final String CAL_PHRED_AVG				=	"avg";
}
