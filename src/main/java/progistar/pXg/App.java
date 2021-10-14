package progistar.pXg;

import progistar.pXg.constants.Parameters;
import progistar.pXg.constants.RunInfo;
import progistar.pXg.data.parser.ParameterParser;
import progistar.pXg.processor.Master;

/**
 * Hello world!
 *
 */
public class App 
{
    public static void main( String[] args )
    {
    	
    	long startTime = System.currentTimeMillis();
    	
    	ParameterParser.parseParams("params.txt");
    	
    	Master.ready(Parameters.genomicAnnotationFilePath, Parameters.sequenceFilePath, Parameters.peptideFilePath);
    	Master.run();
    	
    	long endTime = System.currentTimeMillis();
    	
    	RunInfo.printProcessedChromosomes();
    	System.out.println("\tTotal Elapsed time: "+((endTime-startTime)/1000) + " sec with "+RunInfo.totalProcessedPeptides +" peptides and " +RunInfo.totalProcessedReads+" reads");    	
    }
}
