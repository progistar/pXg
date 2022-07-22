package progistar.pXg;

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
    	// fail to parse params
    	if(ParameterParser.parseParams(args) == -1) {
    		System.exit(1);
    	}
    	
    	Master.ready(Parameters.genomicAnnotationFilePath, Parameters.sequenceFilePath, Parameters.peptideFilePath);
    	Master.run();
    	
    	long endTime = System.currentTimeMillis();
    	
    	RunInfo.printProcessedChromosomes();
    	RunInfo.printPSMCutoff();
    	RunInfo.printFilterStat();
    	RunInfo.printRankPSM();
    	System.out.println("\tTotal Elapsed time: "+((endTime-startTime)/1000) + " sec with "+RunInfo.totalProcessedPeptides +" peptides and " +RunInfo.totalProcessedReads+" reads");
    	
    	Logger.append("\tTotal Elapsed time: "+((endTime-startTime)/1000) + " sec with "+RunInfo.totalProcessedPeptides +" peptides and " +RunInfo.totalProcessedReads+" reads");
    	Logger.newLine();
    	Logger.close();
    }
}