package progistar.thirdparty.sSim;

import java.io.IOException;
import java.util.ArrayList;

public class FindPeakMatch {

	public static void main(String[] args) throws IOException {
		PeptideLoader.loadPeptideList("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/10.Unmodified_10ppm_basic/ProteomeToolsHLA/Matched_NCPeptideList.tsv");
		/*
		 * Load ProteomTools
		 */
		SpectrumGen.loadProteomeTools("/Volumes/One Touch/pxgSelectedRes",
				"/Volumes/One Touch/pxgSelected",
				PeptideLoader.matchedPeptides);
		/*
		SpectrumGen.loadpXg("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/10.Unmodified_10ppm_basic/ProteomeToolsHLA/UniqueNoncanonicalPSMs.tsv",
				"/Users/gistar/projects/pXg/Laumont_NatCommun2016/BLCL_spectra",
				PeptideLoader.matchedPeptides);
		*/
		ArrayList<Spectra> list = SpectrumGen.spectraList;
		System.out.println(list.size()+" were retrived");
		int[] msLevel = {2};
		for(Spectra spectra : list) {
			System.out.println(spectra.getFileName()+": " +spectra.sizeOfSpectra());
			spectra.writeToMGF(spectra.getFileName().replace(".mgf", ".selected.mgf"), msLevel);
		}
	}
}
