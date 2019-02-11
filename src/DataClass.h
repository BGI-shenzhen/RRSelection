
////////////////////////swimming in the sea & flying in the sky //////////////////


/*
 * DataClass.h
 *
 *  Created on: 2011-11-21
 *      Author: hewm@genomics.org.cn
 */

#ifndef DataClass_H_
#define DataClass_H_


using namespace std;


class In3str1v {
	public:
		string InStr1;
		string InStr2;
		string SubPop;
		int  bin ;
		int Masked ;
		int  WindowSize ;
		int  count ;
		bool TF;
		double  Step ;
		double Pvalue;
		In3str1v()
		{
			InStr1="";
			InStr2="";
			SubPop="";
			bin=1;
			Step=0.2;
			count=1;
			WindowSize=300;
			TF=true;
			Masked=10;
			Pvalue=0.005;
		}
};


class Para_18 {
	public:
		string input ;
		string output ;
		double Het ;
		double Miss ;
		double MAF ;
		int Cut3base ;
		Para_18()
		{
			input="";
			output="";
			Het=0.88 ;
			Miss=0.25 ;
			Cut3base=0;
			MAF=0.05;
		}
};


class  StarRsult
{
	public:
		double Count;
		double sumRR;
//		double sumD;
		StarRsult()
		{
			Count=0.0;
			sumRR=0.0;
//			sumD=0.0;
		}
};


class statementVar
{
	public:
		double LN10;
		int Asize ;
		
		double ALL_count ;
		map <string,int>  :: iterator it;
		double D_A;
		double Cal_A ;
		double Cal_B ;
		double D_max;


		double XpA1_pA2;
		double XpA1_pB2;
		double XpB1_pA2;
		double XpB1_pB2;


		double known[5];
		double probHaps[4];
		double lsurface[101];
		double cut5off ;
		unsigned short int DDE[3][3];
		double pA1, pB1, pA2, pB2, loglike1, loglike0;

		double tmpAA, tmpAB, tmpBA, tmpBB, dpr;// tmp2AA, tmp2AB, tmp2BA, tmp2BB;
		int i;
		short int low_i ;
		short int high_i;
		double total_prob;
		double  sum_prob ;
		//		double tmp;//g,h,m,tmp,r;

		statementVar()
		{
			LN10=log(10.0);
			total_prob= 0.0;
			sum_prob=0.0;
			low_i = 0;
			high_i = 0;

		}


};


struct BaseType
{
	unsigned short int Value:2 ;
};


#endif /* DataClass_H_ */

//////////////// swimming in the sky and flying in the sea ////////////////
