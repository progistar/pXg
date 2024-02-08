package progistar.pXg;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import progistar.pXg.constants.Constants;
import progistar.pXg.constants.Parameters;
import progistar.pXg.data.pXgRecord;
import progistar.pXg.data.parser.pXgParser;

public class BuildSequenceDB {

	public static void main(String[] args) throws IOException {
		long startTime = System.currentTimeMillis();
		System.out.println(Constants.VERSION+" "+Constants.RELEASE);
		System.out.println(Constants.INTRODUCE);
		System.out.println();
		
		// parse the options
		parseOptions(args);
		
		ArrayList<pXgRecord> records = pXgParser.parse(new File(Parameters.sequencedbInputPath));
		
		// write the records
		File output = new File(Parameters.sequencedbOutputPath);
		System.out.println("write "+output.getName());
		BufferedWriter BW = new BufferedWriter(new FileWriter(output));
		
		int writeCnt = 0;
		for(int i=0; i<records.size(); i++) {
			pXgRecord record = records.get(i);
			boolean isWrite = false;
			if(Parameters.isIncludedCanonical && record.isCanonical()) {
				isWrite = true;
			}
			else if(Parameters.isIncludedNoncanonical && !record.isCanonical()) {
				isWrite = true;
			}
			
			if(isWrite) {
				writeCnt++;
				BW.append(record.getHeader());
				BW.newLine();
				BW.append(record.getTranslatedSequence());
				BW.newLine();
			}
		}
		System.out.println("A total of "+writeCnt+" were written in the sequence database");
		long endTime = System.currentTimeMillis();
		System.out.println("A total elapsed time: "+(endTime - startTime)/1000 +" sec");
		BW.close();
	}
	
	/**
	 * Parse and apply arguments
	 * 
	 * @param args
	 */
	public static void parseOptions (String[] args) {
		CommandLine cmd = null;
		Options options = new Options();
		
		// Mandatory
		Option optionInput = Option.builder("i")
				.longOpt("input").argName("pXg input")
				.hasArg()
				.required(true)
				.desc("pXg file path")
				.build();
		Option optionOutput = Option.builder("o")
				.longOpt("output").argName("Database output")
				.hasArg()
				.required(true)
				.desc("Sequence database output path")
				.build();
		
		// Optional
		Option optionCanonical = Option.builder("c")
				.longOpt("canonical").argName("Canonical peptides")
				.required(false)
				.desc("Include canonical peptides")
				.build();
		Option optionNoncanonical = Option.builder("n")
				.longOpt("noncanonical").argName("Non-canonical peptides")
				.required(false)
				.desc("Include non-canonical peptides")
				.build();
		Option optionFlank = Option.builder("f")
				.longOpt("flank").argName("Flank sequences")
				.required(false)
				.desc("Include flank sequences")
				.build();
		
		
		options.addOption(optionInput)
		.addOption(optionOutput)
		.addOption(optionCanonical)
		.addOption(optionNoncanonical)
		.addOption(optionFlank);
		
		CommandLineParser parser = new DefaultParser();
	    HelpFormatter helper = new HelpFormatter();
	    boolean isFail = false;
	    
		try {
		    cmd = parser.parse(options, args);
		    
		    if(!cmd.hasOption("c") && !cmd.hasOption("n")) {
		    	isFail = true;
		    	System.out.println("-c or -n must be included.");
		    }
		    if(cmd.hasOption("c")) {
		    	Parameters.isIncludedCanonical = true;
		    	System.out.println("Canonical peptides are included in the sequence database.");
		    }
		    if(cmd.hasOption("n")) {
		    	Parameters.isIncludedNoncanonical = true;
		    	System.out.println("Non-canonical peptides are included in the sequence database.");
		    }
		    if(cmd.hasOption("f")) {
		    	Parameters.isIncludedFlankSequence = true;
		    	System.out.println("Flank sequences are included in the sequence database.");
		    }
		    if(cmd.hasOption("i")) {
		    	Parameters.sequencedbInputPath = cmd.getOptionValue("i");
		    	System.out.println("Input pXg file: "+Parameters.sequencedbInputPath);
		    }
		    if(cmd.hasOption("o")) {
		    	Parameters.sequencedbOutputPath = cmd.getOptionValue("o");
		    	System.out.println("Output the generated sequence database: "+Parameters.sequencedbOutputPath);
		    }
		} catch (ParseException e) {
			System.out.println(e.getMessage());
			isFail = true;
		}
		
		if(isFail) {
		    helper.printHelp("Usage:", options);
		    System.exit(0);
		}
		
		System.out.println();
	}	
}
