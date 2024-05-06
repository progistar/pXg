package progistar.thirdparty.sSim;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import progistar.thirdparty.netMHCpan.NetMHCpanParser;
import progistar.thirdparty.netMHCpan.NetMHCpanResult;

public class SpectralCompToMSP {

	public static void main(String[] args) throws IOException {

		boolean[] selectedSet = {
				false, // B-LCL1
				false, // B-LCL2
				false, // B-LCL3
				false, // B-LCL4
				false, // DOHH2
				false, // HBL1
				false, // SUDHL4
				false, // THP1-1
				false, // THP1-2
				false, // THP1-3

				false, // DI2T
				false, // DI5T
				false, // IN19T
				false, // IN26T
				false, // IN81T
				false, // R1_IN403T
				false, // R2_IN403T
				false, // IN407T
				false, // IN506T
				false, // IN524T
				false, // IN525T
				false, // IN526T
				false, // IN529T
				false, // R1_IN532T
				false, // R2_IN532T
				false, // M004T
				false, // M009T


				true, // DI2N
				true, // DI5N
				true, // IN19N
				true, // IN26N
				true, // IN81N
				true, // IN403N
				true, // IN407N
				true, // IN506N
				true, // IN524N
				true, // IN525N
				true, // IN526N
				true, // IN529N
				true, // IN532N
				true, // M004N
				true // M009N
		};

		String[] exSpectraSet = {
				"/Users/gistar/projects/pXg/MSMS/B_LCL1.mgf",
				"/Users/gistar/projects/pXg/MSMS/B_LCL2.mgf",
				"/Users/gistar/projects/pXg/MSMS/B_LCL3.mgf",
				"/Users/gistar/projects/pXg/MSMS/B_LCL4.mgf",
				"/Users/gistar/projects/pXg/MSMS/DOHH2_400M_050219.mgf",
				"/Users/gistar/projects/pXg/MSMS/HBL1_DMSO_200M_050219.mgf",
				"/Users/gistar/projects/pXg/MSMS/SUDHL4_400M_050219.mgf",
				"/Users/gistar/projects/pXg/MSMS/THP1_1.mgf",
				"/Users/gistar/projects/pXg/MSMS/THP1_2.mgf",
				"/Users/gistar/projects/pXg/MSMS/THP1_3.mgf",

				"/Users/gistar/projects/GastricCancer_NCC/mgf/20221125_KHK_T1_05.scannum.mgf", // DI2T
				"/Users/gistar/projects/GastricCancer_NCC/mgf/20221125_KHK_T1_06.scannum.mgf", // DI5T
				"/Users/gistar/projects/GastricCancer_NCC/mgf/20221125_KHK_T1_12.scannum.mgf", // IN19T
				"/Users/gistar/projects/GastricCancer_NCC/mgf/20221125_KHK_T1_09.scannum.mgf", // IN26T
				"/Users/gistar/projects/GastricCancer_NCC/mgf/20221125_KHK_T1_10.scannum.mgf", // IN81T
				"/Users/gistar/projects/GastricCancer_NCC/mgf/20201110_HLA-I_4_R1.scannum.mgf", // R1_IN403T
				"/Users/gistar/projects/GastricCancer_NCC/mgf/20201110_HLA-I_8_R1.scannum.mgf", // R2_IN403T
				"/Users/gistar/projects/GastricCancer_NCC/mgf/20221125_KHK_T1_11.scannum.mgf", // IN407T
				"/Users/gistar/projects/GastricCancer_NCC/mgf/20221125_KHK_T1_04.scannum.mgf", // IN506T
				"/Users/gistar/projects/GastricCancer_NCC/mgf/20221125_KHK_T1_08.scannum.mgf", // IN524T
				"/Users/gistar/projects/GastricCancer_NCC/mgf/20221125_KHK_T1_01.scannum.mgf", // IN525T
				"/Users/gistar/projects/GastricCancer_NCC/mgf/20201110_HLA-I_3_R1.scannum.mgf", // IN526T
				"/Users/gistar/projects/GastricCancer_NCC/mgf/20221125_KHK_T1_03.scannum.mgf", // IN529T
				"/Users/gistar/projects/GastricCancer_NCC/mgf/20201110_HLA-I_2_R1.scannum.mgf", // R1_IN532T
				"/Users/gistar/projects/GastricCancer_NCC/mgf/20201110_HLA-I_5_R1.scannum.mgf", // R2_IN532T
				"/Users/gistar/projects/GastricCancer_NCC/mgf/20221125_KHK_T1_07.scannum.mgf", // M004T
				"/Users/gistar/projects/GastricCancer_NCC/mgf/230217_MHC_I_160min_M009T.scannum.mgf", // M009T

				"/Users/gistar/projects/GastricCancer_NCC/mgf/20221125_KHK_N1_05.scannum.mgf", // DI2N
				"/Users/gistar/projects/GastricCancer_NCC/mgf/20221125_KHK_N1_06.scannum.mgf", // DI5N
				"/Users/gistar/projects/GastricCancer_NCC/mgf/20221125_KHK_N1_12.scannum.mgf", // IN19N
				"/Users/gistar/projects/GastricCancer_NCC/mgf/20221125_KHK_N1_09.scannum.mgf", // IN26N
				"/Users/gistar/projects/GastricCancer_NCC/mgf/20221125_KHK_N1_10.scannum.mgf", // IN81N
				"/Users/gistar/projects/GastricCancer_NCC/mgf/20201110_HLA-I_1_R1.scannum.mgf", // IN403N
				"/Users/gistar/projects/GastricCancer_NCC/mgf/20221125_KHK_N1_11.scannum.mgf", // IN407N
				"/Users/gistar/projects/GastricCancer_NCC/mgf/20221125_KHK_N1_04.scannum.mgf", // IN506N
				"/Users/gistar/projects/GastricCancer_NCC/mgf/20221125_KHK_N1_08.scannum.mgf", // IN524N
				"/Users/gistar/projects/GastricCancer_NCC/mgf/20221125_KHK_N1_01.scannum.mgf", // IN525N
				"/Users/gistar/projects/GastricCancer_NCC/mgf/20201110_HLA-I_7_R1.scannum.mgf", // IN526N
				"/Users/gistar/projects/GastricCancer_NCC/mgf/20221125_KHK_N1_03.scannum.mgf", // IN529N
				"/Users/gistar/projects/GastricCancer_NCC/mgf/20201110_HLA-I_6_R1.scannum.mgf", // IN532N
				"/Users/gistar/projects/GastricCancer_NCC/mgf/20221125_KHK_N1_07.scannum.mgf", // M004N
				"/Users/gistar/projects/GastricCancer_NCC/mgf/230118_MHC_I_160min_M009N.scannum.mgf" // M009N
		};

		String[] predSpectraSet = {
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/B-LCL1.msp",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/B-LCL2.msp",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/B-LCL3.msp",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/B-LCL4.msp",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/DOHH2.msp",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/HBL1.msp",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/SUDHL4.msp",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/THP1-1.msp",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/THP1-2.msp",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/THP1-3.msp",

				"/Users/gistar/projects/GastricCancer_NCC/pXg/DI2T.msp", // DI2T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/DI5T.msp", // DI5T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN19T.msp", // IN19T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN26T.msp", // IN26T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN81T.msp", // IN81T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/R1_IN403T.msp", // R1_IN403T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/R2_IN403T.msp", // R2_IN403T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN407T.msp", // IN407T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN506T.msp", // IN506T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN524T.msp", // IN524T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN525T.msp", // IN525T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN526T.msp", // IN526T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN529T.msp", // IN529T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/R1_IN532T.msp", // R1_IN532T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/R2_IN532T.msp", // R2_IN532T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/M004T.msp", // M004T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/M009T.msp", // M009T

				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/DI2N.msp", // DI2N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/DI5N.msp", // DI5N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN19N.msp", // IN19N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN26N.msp", // IN26N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN81N.msp", // IN81N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN403N.msp", // IN403N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN407N.msp", // IN407N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN506N.msp", // IN506N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN524N.msp", // IN524N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN525N.msp", // IN525N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN526N.msp", // IN526N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN529N.msp", // IN529N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN532N.msp", // IN532N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/M004N.msp", // M004N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/M009N.msp" // M009N
		};

		String[] deepLCSet = {
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/B-LCL1.deeplc",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/B-LCL2.deeplc",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/B-LCL3.deeplc",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/B-LCL4.deeplc",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/DOHH2.deeplc",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/HBL1.deeplc",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/SUDHL4.deeplc",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/THP1-1.deeplc",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/THP1-2.deeplc",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/THP1-3.deeplc",

				"/Users/gistar/projects/GastricCancer_NCC/pXg/DI2T.deeplc", // DI2T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/DI5T.deeplc", // DI5T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN19T.deeplc", // IN19T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN26T.deeplc", // IN26T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN81T.deeplc", // IN81T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/R1_IN403T.deeplc", // R1_IN403T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/R2_IN403T.deeplc", // R2_IN403T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN407T.deeplc", // IN407T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN506T.deeplc", // IN506T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN524T.deeplc", // IN524T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN525T.deeplc", // IN525T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN526T.deeplc", // IN526T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN529T.deeplc", // IN529T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/R1_IN532T.deeplc", // R1_IN532T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/R2_IN532T.deeplc", // R2_IN532T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/M004T.deeplc", // M004T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/M009T.deeplc", // M009T

				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/DI2N.deeplc", // DI2N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/DI5N.deeplc", // DI5N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN19N.deeplc", // IN19N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN26N.deeplc", // IN26N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN81N.deeplc", // IN81N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN403N.deeplc", // IN403N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN407N.deeplc", // IN407N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN506N.deeplc", // IN506N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN524N.deeplc", // IN524N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN525N.deeplc", // IN525N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN526N.deeplc", // IN526N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN529N.deeplc", // IN529N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN532N.deeplc", // IN532N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/M004N.deeplc", // M004N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/M009N.deeplc" // M009N

		};

		String[] netmhcpanSet = {
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/B-LCL1.netmhcpan.xls",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/B-LCL2.netmhcpan.xls",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/B-LCL3.netmhcpan.xls",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/B-LCL4.netmhcpan.xls",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/DOHH2.netmhcpan.xls",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/HBL1.netmhcpan.xls",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/SUDHL4.netmhcpan.xls",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/THP1-1.netmhcpan.xls",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/THP1-2.netmhcpan.xls",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/THP1-3.netmhcpan.xls",

				"/Users/gistar/projects/GastricCancer_NCC/pXg/DI2T.netmhcpan.xls", // DI2T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/DI5T.netmhcpan.xls", // DI5T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN19T.netmhcpan.xls", // IN19T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN26T.netmhcpan.xls", // IN26T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN81T.netmhcpan.xls", // IN81T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/R1_IN403T.netmhcpan.xls", // R1_IN403T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/R2_IN403T.netmhcpan.xls", // R2_IN403T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN407T.netmhcpan.xls", // IN407T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN506T.netmhcpan.xls", // IN506T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN524T.netmhcpan.xls", // IN524T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN525T.netmhcpan.xls", // IN525T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN526T.netmhcpan.xls", // IN526T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN529T.netmhcpan.xls", // IN529T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/R1_IN532T.netmhcpan.xls", // R1_IN532T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/R2_IN532T.netmhcpan.xls", // R2_IN532T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/M004T.netmhcpan.xls", // M004T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/M009T.netmhcpan.xls", // M009T

				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/DI2N.netmhcpan.xls", // DI2N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/DI5N.netmhcpan.xls", // DI5N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN19N.netmhcpan.xls", // IN19N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN26N.netmhcpan.xls", // IN26N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN81N.netmhcpan.xls", // IN81N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN403N.netmhcpan.xls", // IN403N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN407N.netmhcpan.xls", // IN407N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN506N.netmhcpan.xls", // IN506N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN524N.netmhcpan.xls", // IN524N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN525N.netmhcpan.xls", // IN525N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN526N.netmhcpan.xls", // IN526N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN529N.netmhcpan.xls", // IN529N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN532N.netmhcpan.xls", // IN532N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/M004N.netmhcpan.xls", // M004N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/M009N.netmhcpan.xls" // M009N
		};

		String[] pinSet = {
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/B-LCL1.pxg.pin",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/B-LCL2.pxg.pin",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/B-LCL3.pxg.pin",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/B-LCL4.pxg.pin",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/DOHH2.pxg.pin",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/HBL1.pxg.pin",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/SUDHL4.pxg.pin",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/THP1-1.pxg.pin",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/THP1-2.pxg.pin",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/THP1-3.pxg.pin",

				"/Users/gistar/projects/GastricCancer_NCC/pXg/DI2T.pxg.pin", // DI2T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/DI5T.pxg.pin", // DI5T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN19T.pxg.pin", // IN19T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN26T.pxg.pin", // IN26T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN81T.pxg.pin", // IN81T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/R1_IN403T.pxg.pin", // R1_IN403T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/R2_IN403T.pxg.pin", // R2_IN403T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN407T.pxg.pin", // IN407T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN506T.pxg.pin", // IN506T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN524T.pxg.pin", // IN524T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN525T.pxg.pin", // IN525T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN526T.pxg.pin", // IN526T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN529T.pxg.pin", // IN529T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/R1_IN532T.pxg.pin", // R1_IN532T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/R2_IN532T.pxg.pin", // R2_IN532T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/M004T.pxg.pin", // M004T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/M009T.pxg.pin", // M009T

				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/DI2N.pxg.pin", // DI2N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/DI5N.pxg.pin", // DI5N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN19N.pxg.pin", // IN19N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN26N.pxg.pin", // IN26N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN81N.pxg.pin", // IN81N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN403N.pxg.pin", // IN403N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN407N.pxg.pin", // IN407N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN506N.pxg.pin", // IN506N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN524N.pxg.pin", // IN524N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN525N.pxg.pin", // IN525N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN526N.pxg.pin", // IN526N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN529N.pxg.pin", // IN529N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN532N.pxg.pin", // IN532N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/M004N.pxg.pin", // M004N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/M009N.pxg.pin" // M009N

		};

		String[] pxgSet = {
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/B-LCL1.pxg",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/B-LCL2.pxg",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/B-LCL3.pxg",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/B-LCL4.pxg",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/DOHH2.pxg",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/HBL1.pxg",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/SUDHL4.pxg",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/THP1-1.pxg",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/THP1-2.pxg",
				"/Users/gistar/eclipse-workspace/pXg/test/high_score_decoy/nocut/THP1-3.pxg",

				"/Users/gistar/projects/GastricCancer_NCC/pXg/DI2T.pxg", // DI2T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/DI5T.pxg", // DI5T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN19T.pxg", // IN19T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN26T.pxg", // IN26T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN81T.pxg", // IN81T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/R1_IN403T.pxg", // R1_IN403T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/R2_IN403T.pxg", // R2_IN403T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN407T.pxg", // IN407T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN506T.pxg", // IN506T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN524T.pxg", // IN524T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN525T.pxg", // IN525T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN526T.pxg", // IN526T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/IN529T.pxg", // IN529T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/R1_IN532T.pxg", // R1_IN532T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/R2_IN532T.pxg", // R2_IN532T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/M004T.pxg", // M004T
				"/Users/gistar/projects/GastricCancer_NCC/pXg/M009T.pxg", // M009T

				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/DI2N.pxg", // DI2N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/DI5N.pxg", // DI5N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN19N.pxg", // IN19N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN26N.pxg", // IN26N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN81N.pxg", // IN81N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN403N.pxg", // IN403N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN407N.pxg", // IN407N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN506N.pxg", // IN506N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN524N.pxg", // IN524N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN525N.pxg", // IN525N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN526N.pxg", // IN526N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN529N.pxg", // IN529N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/IN532N.pxg", // IN532N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/M004N.pxg", // M004N
				"/Users/gistar/projects/GastricCancer_NCC_2nd_analysis/pXg_normal/M009N.pxg" // M009N

		};


		for(int idx = 0; idx<pinSet.length; idx++) {
			if(!selectedSet[idx]) {
				continue;
			}

			System.out.println(pinSet[idx]);
			System.out.println(exSpectraSet[idx]);
			System.out.println(predSpectraSet[idx]);
			System.out.println(deepLCSet[idx]);
			System.out.println(netmhcpanSet[idx]);
			// index by title
			Spectra exSpectra = new Spectra(exSpectraSet[idx], Spectra.FILE_TYPE_MGF);
			// index by peptide/charge
			Spectra predSpectra = new Spectra(predSpectraSet[idx], Spectra.FILE_TYPE_MSP);
			// deepLC best
			File deepLCRes = new File(deepLCSet[idx]);
			// HLA binding prediction
			File netmhcpan = new File(netmhcpanSet[idx]);

			NetMHCpanResult netMHCpanResult = NetMHCpanParser.parseNetMHCpan(netmhcpan.getAbsolutePath());

			File pinFile = new File(pinSet[idx]);
			String line = null;

			BufferedReader BRlc = new BufferedReader(new FileReader(deepLCRes));
			Hashtable<String, Double> bestDelta = new Hashtable<>();
			Hashtable<String, ArrayList<Double>> deltaRTs = new Hashtable<>();
			BRlc.readLine();// skip header

			while((line = BRlc.readLine()) != null) {
				String[] fields = line.split("\\,");
				String peptide =fields[1];
				// if observed is "-"
				if(fields[3].equalsIgnoreCase("-")) {
					fields[3] = fields[4];
				}
				Double obRT = Double.parseDouble(fields[3]);
				Double predRT = Double.parseDouble(fields[4]);
				double delta = obRT - predRT;
				Double thisDelta = bestDelta.get(peptide);
				ArrayList<Double> deltaList = deltaRTs.get(peptide);
				if(thisDelta == null || Math.abs(thisDelta) > Math.abs(delta)) {
					bestDelta.put(peptide, delta);
				}

				if(deltaList == null) {
					deltaList = new ArrayList<>();
				}
				deltaList.add(delta);
				deltaRTs.put(peptide, deltaList);
			}

			BRlc.close();

			BufferedReader BR = new BufferedReader(new FileReader(pxgSet[idx]));
			ArrayList<String> pxgRecords = new ArrayList<>();
			String pXgHeader = BR.readLine();
			while((line = BR.readLine()) != null) {
				pxgRecords.add(line);
			}

			BR.close();

			BR = new BufferedReader(new FileReader(pinFile));
			BufferedWriter BW = new BufferedWriter(new FileWriter(pinFile.getAbsolutePath().replace(".pin", ".predfeat.pin")));
			BufferedWriter BWpxg = new BufferedWriter(new FileWriter(pxgSet[idx].replace(".pxg", ".predfeat.pxg")));
			BWpxg.append(pXgHeader+"\tSA\tBestDeltaRT\tmLog2BestELRank");
			BWpxg.newLine();


			StringBuilder output = new StringBuilder();
			String[] headers = BR.readLine().split("\t");
			for(int i=0; i<headers.length-2; i++) {
				output.append(headers[i]).append("\t");
			}
			output.append("SA\tBestDeltaRT\tmLog2BestELRank\t").append(headers[headers.length-2]).append("\t").append(headers[headers.length-1]);
			BW.append(output.toString());
			BW.newLine();
			int pxgIdx = 0;
			while((line = BR.readLine()) != null) {
				output.setLength(0);
				String[] fields = line.split("\t");
				String uniqueID = fields[0];
				String fileName = uniqueID.split("\\|")[0].replace(".scannum.mgf", "");
				String scanNum = uniqueID.split("\\|")[1];
				String charge = uniqueID.split("\\|")[2];
				String title = fileName+"."+scanNum+"."+scanNum+"."+charge;
				String peptide = fields[fields.length-2];

				Spectrum s1 = exSpectra.getSpectrumByScanNum(title);
				Spectrum s2 = predSpectra.getSpectrumByScanNum(peptide+"/"+charge);

				s1.setPeptide(peptide);
				s2.setPeptide(peptide);

				double scaScore = SpectralScores.spectralContrastAngle(s1, s2, 0.02, false);

				for(int i=0; i<fields.length-2; i++) {
					output.append(fields[i]).append("\t");
				}

				ArrayList<Double> deltaList = deltaRTs.get(peptide);
				double avgDelta = 0;
				for(Double delta : deltaList) {
					avgDelta += delta;
				}
				avgDelta /= deltaList.size();

				double elrank = -Math.log(netMHCpanResult.peptideToRecord.get(peptide).getBestScore()/2) / Math.log(2);

				output.append(scaScore).append("\t")
				.append(Math.abs(bestDelta.get(peptide))).append("\t")
				.append(elrank).append("\t")
				.append(fields[fields.length-2]).append("\t").append(fields[fields.length-1]);
				BW.append(output.toString());
				BW.newLine();

				BWpxg.append(pxgRecords.get(pxgIdx++)).append("\t").append(scaScore+"").append("\t")
				.append(Math.abs(bestDelta.get(peptide))+"").append("\t")
				.append(elrank+"");
				BWpxg.newLine();

			}
			BW.close();
			BR.close();
			BWpxg.close();
		}

	}


}
