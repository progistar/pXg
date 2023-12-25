package progistar.pXg;

import progistar.pXg.constants.Constants;
import progistar.pXg.constants.Parameters;
import progistar.pXg.constants.RunInfo;
import progistar.pXg.data.parser.ParameterParser;
import progistar.pXg.processor.Master;
import progistar.pXg.utils.Logger;

/**
 * Hello world!
 *
 */
public class Starter 
{
    public static void main( String[] args )
    {
    	
    	long startTime = System.currentTimeMillis();
    	// set this application id by start time.
    	// the start time is supposed to be an unique id.
    	// this unique id is used to name the temporary files.
    	Constants.UNIQUE_RUN_ID = startTime;
    	
    	// fail to parse params
    	if(ParameterParser.parseParams(args) == -1) {
    		System.exit(1);
    	}
    	
    	Master.ready(Parameters.genomicAnnotationFilePath, Parameters.sequenceFilePath, Parameters.peptideFilePath);
    	Master.run();
    	
    	long endTime = System.currentTimeMillis();
    	
    	RunInfo.printProcessedChromosomes();
    	RunInfo.printFilterStat();
    	System.out.println("\tTotal Elapsed time: "+((endTime-startTime)/1000) + " sec with "+RunInfo.totalProcessedPeptides +" peptides and " +RunInfo.totalProcessedReads+" reads");
    	
    	Logger.append("\tTotal Elapsed time: "+((endTime-startTime)/1000) + " sec with "+RunInfo.totalProcessedPeptides +" peptides and " +RunInfo.totalProcessedReads+" reads");
    	Logger.newLine();
    	Logger.close();
    	
    	// TODO: third-party runner
    	
    }
}