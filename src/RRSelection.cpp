#ifndef RRSelection_H_
#define RRSelection_H_

#include "HeadIN.h"
#include "FileDeal.h"
#include "Calculate.h"
#include "CalOnePop.h"
#include "CalTwoPop.h"

using namespace std;

void  RRSelection_help()
{
	cout <<""
		"\n"
		"\tUsage: RRSelection  -InVCF  <in.vcf.gz>  -OutPut <outPrefix>\n"
		"\n"
		"\t\t-InVCF      <str>     Input SNP VCF Format\n"
		"\t\t-OutPut     <str>     OutPut sliding stat mean r^2 Result\n"
		"\n"
		"\t\t-SubGroup   <str>     one/two sub-group Sample List File,-h for more help\n"
		"\t\t-Windows    <int>     Sliding windows bin (kb),MaxDis between two pairwise SNP[300]\n"
		"\t\t-Step       <float>   Step ratio(0,1] of windows,1:NoOverlap [0.2]\n"
		"\t\t-Masked     <int>     Masked windows when the SNP Number too low[10]\n"
		"\t\t\n"
		"\t\t-MAF        <float>   Min minor allele frequency filter [0.05]\n"
		"\t\t-Het        <float>   Max ratio of het allele filter [0.88]\n"
		"\t\t-Miss       <float>   Max ratio of miss allele filter [0.25]\n"
		"\t\t\n"
		"\t\t-Pvalue     <float>   T-test Pvalue to pick out selection region[0.005]\n"
		"\t\t-KeepR                Keep Rscript used to modify and plots\n"
		"\t\t\n"
		"\t\t-help                 See more help [hewm2008 Beta v0.85]\n"
		"\n";
}

void More_HelpRRSelection()
{
	cout<<""
		"\n"
		"\t\t More Help document please see the Manual.pdf file\n"
		"\t\t 1. Para [-SubGroup] only can accept one or two groups\n"
		"\t\t   (1.1)  For two sub group, the input file should be fixed format(groupID : samaple list), For example :\n"
		"\t\t            cultivation :  CulSampleNameA CulSampleNameB ... CulSampleNameN\n"
		"\t\t            wild :  WildSampleNameA WildSampleNameB ... WildSampleNameN\n"
		"\t\t   (1.2)  For one sub group, the input file should be fixed or not. Rows or Columns is ok\n"
		"\t\t            wild :  WildSampleNameA WildSampleNameB ... WildSampleNameN\n"
		"\n";
}

int LDdecay_help01(int argc, char **argv , In3str1v * paraFA04, Para_18 * para_18)
{
	if (argc <2 ) {RRSelection_help();return 0;}

	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InVCF" ||  flag  == "i")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr1=argv[i];
		}
		else if (flag  == "SubGroup"  ||  flag  == "s")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->SubPop=argv[i];
		}
		else if (flag  ==  "OutPut" ||  flag  == "o" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr2=argv[i];
		}
		else if (flag  == "Het" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			para_18->Het=atof(argv[i]);
		}
		else if (flag  == "Masked")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->Masked=atoi(argv[i]);
		}
		else if (flag == "MAF")
		{
			if(i + 1== argc) {LogLackArg(flag);return 0;}
			i++;
			para_18->MAF=atof(argv[i]);
		}
		else if (flag == "Pvalue")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			paraFA04->Pvalue=atof(argv[i]);
		}
		else if (flag == "Miss")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_18->Miss=atof(argv[i]);
		}
		else if (flag == "Windows" )
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			paraFA04->WindowSize=atoi(argv[i]);
		}
		else if (flag == "Step")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			paraFA04->Step=atof(argv[i]);
			if (  (paraFA04->Step)>0  &&  (paraFA04->Step)<=1.0) {} else
			{
				cerr  <<"ratio should be in the Number region (0,1]\n";
				return 0;
			}
		}
		else if (flag == "KeepR")
		{
			paraFA04->TF=false;
		}
		else if (flag == "help"  ||  flag == "h")
		{
			More_HelpRRSelection();return 0;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}

	if  ( (paraFA04->InStr2).empty() ||  (paraFA04->InStr1).empty()  )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}


	string Stat=(paraFA04->InStr2);
	string ext =Stat.substr(Stat.rfind('.') ==string::npos ? Stat.length() : Stat.rfind('.') + 1);
	if (ext == "gz")
	{
		(paraFA04->InStr2)=(paraFA04->InStr2).substr(0,(paraFA04->InStr2).length()-3);
	}

	Stat=(paraFA04->InStr2);
	ext =Stat.substr(Stat.rfind('/') ==string::npos ? Stat.length() : Stat.rfind('/') + 1);

	if (ext != "stat")
	{
		ext =Stat.substr(Stat.rfind('.') ==string::npos ? Stat.length() : Stat.rfind('.') + 1);
		if (ext == "stat")
		{
			(paraFA04->InStr2)=(paraFA04->InStr2).substr(0,(paraFA04->InStr2).length()-5);
		}
	}
	return 1 ;
}



int main(int argc, char *argv[])
{
	In3str1v *paraFA04 = new In3str1v;
	Para_18 * para_18 = new Para_18 ;
	if ((LDdecay_help01(argc, argv, paraFA04, para_18)==0))
	{
		delete paraFA04 ;
		delete para_18 ;
		return 1 ;
	}

	(paraFA04->WindowSize)=(paraFA04->WindowSize)*1000;
	(paraFA04->count)=int(1/(paraFA04->Step));
	(paraFA04->Step)=1.0/(paraFA04->count);
	(paraFA04->bin)=int((paraFA04->WindowSize)/(paraFA04->count));
	(paraFA04->Masked)=int(((paraFA04->Masked)*((paraFA04->Masked)-1)/2)* 1.6449); //  sum(1/(n*n))= pi^2/6 =1.6449

	//*///////////////////////////Test  Out File is OK ///////////////////////////////*//
	string Stat=(paraFA04->InStr2);
	Stat=(paraFA04->InStr2)+".winRR.gz";
	ogzstream OUTTest ((Stat).c_str());
	if((!OUTTest.good()))
	{
		cerr << "open OUT File error: "<<Stat<<endl;
		delete para_18;
		delete  paraFA04 ; return 1;
	}
	OUTTest.close();


	map <string,map <llong, vector <BaseType> > >  SNPList ;
	map <string,llong> MaxSNPSite;
	map <string,map <llong, vector <BaseType> > >  SNPListSub ;
	int Flag_for_pro=0;

	map <string,int >  GroupID ;

	/////*   VCF   IN Deal //////////////////////*////
	bool PhaseVCF=DeleVCFPhase(paraFA04->InStr1);

	int IFTwoGroup=1 ;
	if (!(paraFA04->SubPop).empty() )
	{
		map <string,int >  SubVetor ;
		IFTwoGroup=Identify_subgroup(paraFA04->SubPop, SubVetor,GroupID );
		if ( IFTwoGroup ==2 )   // two group Read
		{
			if   ( PhaseVCF )
			{
				Read_TwoPopVCF_IN_Phase( paraFA04, para_18 ,SubVetor , SNPList,SNPListSub ,MaxSNPSite, Flag_for_pro);
			}
			else
			{
				Read_TwoPopVCF_IN( paraFA04, para_18 ,SubVetor, SNPList,SNPListSub,MaxSNPSite, Flag_for_pro);
			}
		}
		else if  (IFTwoGroup  >2 ) // muti group  wrong  exit
		{
			delete para_18;
			delete paraFA04 ; return 1 ;			
		}
		else  // one group with sub pop
		{
			if  ( PhaseVCF )
			{
				Read_SubPopVCF_IN_Phase( paraFA04, para_18 , SNPList,MaxSNPSite, Flag_for_pro);
			}
			else
			{
				Read_SubPopVCF_IN( paraFA04, para_18 , SNPList,MaxSNPSite, Flag_for_pro);
			}
		}
	}
	else // one group with all sample
	{
		if  ( PhaseVCF )
		{
			Read_VCF_IN_Phase( paraFA04, para_18 , SNPList, MaxSNPSite, Flag_for_pro);
		}
		else
		{
			Read_VCF_IN( paraFA04, para_18 , SNPList,MaxSNPSite, Flag_for_pro);
		}
	}



	//*///////////////////////////PairWise Compare//////////////////////////////////*//

	if ( IFTwoGroup ==1 ) // run one group 
	{
		cout<<"##begin pair-wise R^2 cal. after filter Remain SNP Number :\t"<<Flag_for_pro<<"\n#\% windows bin is\t"<<(paraFA04->WindowSize)<<endl;
		SilidingRRCal_One(paraFA04, para_18 , SNPList , MaxSNPSite );
	}
	else if ( IFTwoGroup ==2 )// run two group
	{
		cout<<"##begin pair-wise R^2 cal. after filter Remain SNP Number :\t"<<Flag_for_pro<<"\n#\% windows bin is\t"<<(paraFA04->WindowSize)<<endl;
		SilidingRRCal_Tow(paraFA04, para_18 , SNPList , SNPListSub , MaxSNPSite , GroupID );
	}

	delete para_18 ;
	delete paraFA04 ;
	return 0;

}


#endif // RRSelection_H_ //
///////// swimming in the sky and flying in the sea ////////////


