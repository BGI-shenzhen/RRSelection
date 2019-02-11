#ifndef processOnePop_H_
#define processOnePop_H_

#include "Bmath/pnorm.c"

using namespace std;

//////////////////////////////////////////////////////////////////  one Group  Deal ///////////////////////////////////////////////////


void DrawRRSliding_One( string File , In3str1v * paraFA04, Para_18 * para_18 )
{
	//  OUTSS<<"Chr\tStart\tEnd\tRLV\tSum_r^2\tCount\tZScore\tPvalue\n";

	string OutPlotTmp=(paraFA04->InStr2)+".tmp";
	ofstream  OUTPLOT (OutPlotTmp.c_str());

	igzstream LIST (File.c_str(),ifstream::in);

	string chr="";
	getline(LIST,chr);
	getline(LIST,chr);
	string Start,End,RLV;
	map <string , int >  VecChr;
	map <string , map <int , string > >  Value;
	map <string , map <int , string > > :: iterator iTValue;

	map  <string , int >  :: iterator itVecChr ;
	while(!LIST.eof())
	{
		string  line ;
		getline(LIST,line);
		if (line.length()<=0) { continue ;}
		istringstream isone (line,istringstream::in);
		isone>>chr>>Start>>End>>RLV ;
		itVecChr=VecChr.find(chr);
		if (itVecChr==VecChr.end())
		{
			VecChr.insert(map <string,int>::value_type(chr,1));
			map <int , string >  DataEE;
			DataEE[0]=RLV;
			Value.insert(map <string , map <int , string > > ::value_type(chr,DataEE));
		}
		else
		{
			iTValue=Value.find(chr);
			(iTValue->second).insert(map <int , string >::value_type(itVecChr->second,RLV));
			itVecChr->second++;
		}
	}
	LIST.close();

	int MaxXX=0;
	OUTPLOT<<"Site";
	for(itVecChr=VecChr.begin();itVecChr!=VecChr.end();itVecChr++)
	{
		if (MaxXX<(itVecChr->second)) {MaxXX=itVecChr->second;}
		OUTPLOT<<"\t"<<itVecChr->first;
	}
	OUTPLOT<<endl;


	map <int , string > ::  iterator  itEEEE;
	for (int ii=0 ; ii< MaxXX ; ii++)
	{
		OUTPLOT<<ii;
		for(itVecChr=VecChr.begin();itVecChr!=VecChr.end();itVecChr++)
		{

			if (ii> (itVecChr->second))
			{
				OUTPLOT<<"\t0.0";
				continue ;
			}
			iTValue=Value.find(itVecChr->first);
			itEEEE=(iTValue->second).find(ii);
			if (itEEEE==(iTValue->second).end())
			{
				OUTPLOT<<"\tNA";
			}
			else
			{
				OUTPLOT<<"\t"<<itEEEE->second ;
			}
		}
		OUTPLOT<<endl;
	}


	OUTPLOT.close();

	string OutPlotr=(paraFA04->InStr2)+".tmp.r";
	ofstream  OUTR (OutPlotr.c_str());

	OUTR<<""
		"\n"
		"library(ggplot2);\nlibrary(reshape2);\n"
		"read.table(\""<<OutPlotTmp<<"\",header=T)->data;\n"
		"data2<-melt(data,id=1);\n"
		"p <- ggplot(data = data2)+geom_tile(aes(x = Site,y = variable, fill=value))+scale_fill_gradient(limit=c(0,1),low ='white',high ='red',na.value=\"grey50\")\n"
		"pdf(\""<<(paraFA04->InStr2)<<".sliding.pdf\",h=6,w=25);\n"
		"p\n"
		"dev.off()\n"
		"png(\""<<(paraFA04->InStr2)<<".sliding.png\",h=6,w=25);\n"
		"p\n"
		"dev.off()\n"<<endl;
	OUTR.close();

	char   buf[2048]={'\0'};
	string cc="which  Rscript  2> /dev/null ";
	memset( buf, '\0', sizeof(buf) );
	FILE   *stream ;
	stream=popen(cc.c_str(),"r") ;
	fread( buf, sizeof(char), sizeof(buf), stream);
	string binPath=buf;
	binPath=binPath.substr(0,binPath.length()-1);
	if (binPath == "" )
	{
		cout <<"\twarning: can't find the [Rscript] in your $PATH ; no png Figure Out"<<endl;
		cout <<"\t\tRscript "<<OutPlotr<<endl;
	}
	else
	{
		if (paraFA04->TF)
		{
			//cc=binPath+"\t"+OutPlotr ;
			cc=binPath+"\t"+OutPlotr+"  ; rm  -rf  "+ OutPlotr+" "+ OutPlotTmp ;
			std::system(cc.c_str()) ;
		}
		else
		{
			cc=binPath+"\t"+OutPlotr ;
			std::system(cc.c_str()) ;
			cout <<"\t\tRePlot : Rscript  "<<OutPlotr<<endl;
		}
	}
}





void DrawRRDis_One( string File1,  string File2, In3str1v * paraFA04, Para_18 * para_18)
{

	string OutPlotTmp=(paraFA04->InStr2)+".tmp2";
	ofstream  OUTPLOT (OutPlotTmp.c_str());

	igzstream LIST (File1.c_str(),ifstream::in);

	string chr="";
	getline(LIST,chr);
	getline(LIST,chr);
	string Start,End,RLVV1;

	map < string,pair <int,int> > ValueV1;
	map < string,pair <int,int> > :: iterator AA_it;

	while(!LIST.eof())
	{
		string  line ;
		getline(LIST,line);
		if (line.length()<=0) {continue;}
		istringstream isone (line,istringstream::in);
		isone>>chr>>Start>>End>>RLVV1;
		if(RLVV1!="NA")
		{
			RLVV1=RLVV1.substr(0,4);
			AA_it=ValueV1.find(RLVV1);
			if (AA_it==ValueV1.end())
			{
				pair <int,int> dataS (1,1);
				ValueV1.insert( map < string,pair <int,int>  >::value_type(RLVV1,dataS));
			}
			else
			{				
				(AA_it->second).first++;
			}
		}
	}
	LIST.close();

	igzstream LIST2 (File2.c_str(),ifstream::in);
	getline(LIST2,chr);
	while(!LIST2.eof())
	{
		string  line ;
		getline(LIST2,line);
		if (line.length()<=0) {continue;}
		istringstream isone (line,istringstream::in);
		isone>>chr>>Start>>End>>RLVV1;
		if(RLVV1!="NA" )
		{
			RLVV1=RLVV1.substr(0,4);
			AA_it=ValueV1.find(RLVV1);
			if(RLVV1<"0.200")
			{
				(AA_it->second).second=2;
			}
			else
			{		
				(AA_it->second).second=3;
			}
		}
	}
	LIST2.close();


	int MaxXX=0;
	OUTPLOT<<"RLV\tNum\tType\n";

	for(AA_it=ValueV1.begin();AA_it!=ValueV1.end();AA_it++)
	{

		if((AA_it->second).second==1)
		{
			OUTPLOT<<AA_it->first<<"\t"<<(AA_it->second).first<<"\tNoSele\n";
		}
		else if ((AA_it->second).second==2)
		{
			OUTPLOT<<AA_it->first<<"\t"<<(AA_it->second).first<<"\tSeleLow\n";
		}
		else
		{
			OUTPLOT<<AA_it->first<<"\t"<<(AA_it->second).first<<"\tSeleHigh\n";
		}
	}
	OUTPLOT.close();

	string OutPlotr=(paraFA04->InStr2)+".tmp2.r";
	ofstream  OUTR (OutPlotr.c_str());

	OUTR<<""
		"\n"
		"library(ggplot2);library(ggthemes);\n"
		"read.table(\""<<OutPlotTmp<<"\",header=T)->data;\n"
		"RLV<-data[,1];\nNum<-data[,2];\nsele<-data[,3];\n"
		"df <- data.frame(RLV,Num);\n"
		"p <- ggplot(df, aes(x = RLV, y = Num, col = sele,fill=factor(sele) ) ) + geom_bar(stat = \"identity\") + theme_economist();\n"
		"pdf(\""<<(paraFA04->InStr2)<<".Dis.pdf\");\n"
		"p\n"
		"dev.off()\n"
		"png(\""<<(paraFA04->InStr2)<<".Dis.png\");\n"
		"p\n"
		"dev.off()\n"<<endl;
	OUTR.close();




	char   buf[2048]={'\0'};
	string cc="which  Rscript  2> /dev/null ";
	memset( buf, '\0', sizeof(buf) );
	FILE   *stream ;
	stream=popen(cc.c_str(),"r") ;
	fread( buf, sizeof(char), sizeof(buf), stream);
	string binPath=buf;
	binPath=binPath.substr(0,binPath.length()-1);
	if (binPath == "" )
	{
		cout <<"\twarning: can't find the [Rscript] in your $PATH ; no png Figure Out"<<endl;
		cout <<"\t\tRscript "<<OutPlotr<<endl;
	}
	else
	{
		if (paraFA04->TF)
		{
			//cc=binPath+"\t"+OutPlotr ;
			cc=binPath+"\t"+OutPlotr+"  ; rm  -rf  "+ OutPlotr+ " "+ OutPlotTmp;
			std::system(cc.c_str()) ;
		}
		else
		{
			cc=binPath+"\t"+OutPlotr ;
			std::system(cc.c_str()) ;
			cout <<"\t\tRePlot : Rscript  "<<OutPlotr<<endl;
		}
	}
}






int SilidingRRCal_One( In3str1v * paraFA04, Para_18 * para_18 , map <string,map <llong, vector <BaseType> > >  & SNPList ,  map <string,llong>  & MaxSNPSite)
{
	string Stat=(paraFA04->InStr2)+".winRR.temp.gz";
	ogzstream OUT ((Stat).c_str());
	if((!OUT.good()))
	{
		cerr << "open OUT File error: "<<Stat<<endl;
		return 0;
	}
	OUT<<"#Chr\tStart\tEnd\tRLV\tSum_r^2\tCount\n";
	//#####################

	map <string, llong> :: iterator  MaxIT ;
	map <string,map<llong, vector <BaseType>  > > :: iterator key1_it ;

	int  ALL_count=0;
	double SumRLV=0.0;

	for (MaxIT=MaxSNPSite.begin(); MaxIT!=MaxSNPSite.end();MaxIT++)
	{
		int MaxBin=int((MaxIT->second)/(paraFA04->bin))+1;
		(paraFA04->count)++;
		StarRsult **All_Stat = new StarRsult *[(paraFA04->count)];
		for (int kk=0 ; kk<(paraFA04->count) ; kk++)
		{
			All_Stat[kk]=new StarRsult[MaxBin];
		}

		key1_it=SNPList.find(MaxIT->first);
		(paraFA04->count)--;
		PairWiseRRCal( paraFA04,  para_18, key1_it->second ,  All_Stat);

		MaxBin=MaxBin-(paraFA04->count)+1;
		for (int ii=0 ; ii<MaxBin ; ii++)
		{
			double count=0;
			double SumRR=0;
			for (int jj=0; jj<(paraFA04->count); jj++)
			{
				for (int uu=jj  ; uu<(paraFA04->count) ; uu++)
				{
					count+=(All_Stat[uu][ii+jj]).Count;
					SumRR+=(All_Stat[uu][ii+jj]).sumRR;
				}
			}
			if (count<(paraFA04->Masked))  {continue;}
			double RLV=SumRR/count;
			llong Start=ii*(paraFA04->bin);
			llong End=Start+(paraFA04->WindowSize);
			SumRLV+=RLV;
			ALL_count++;
			OUT<<MaxIT->first<<"\t"<<Start<<"\t"<<End<<"\t"<<setprecision(4)<<setiosflags(ios::left)<<setiosflags(ios::fixed)<<RLV<<"\t"<<SumRR<<"\t"<<count<<"\n";
		}

		for (int kk=0 ; kk<(paraFA04->count) ; kk++)
		{
			delete [] All_Stat[kk];
		}

		delete [] All_Stat ;
	}
	OUT.close();



	double meaWolRR=SumRLV/ALL_count;


	igzstream INAA (Stat.c_str(),ifstream::in);
	if (INAA.fail())
	{
		cerr << "open Sub Group IN File error: "<<Stat<<endl;
		return  0;
	}


	string StatSS=(paraFA04->InStr2)+".winRR.gz";
	ogzstream OUTSS ((StatSS).c_str());
	if((!OUTSS.good()))
	{
		cerr << "open OUT File error: "<<Stat<<endl;
		return 0;
	}

	string StatSELE=(paraFA04->InStr2)+".selection.gz";
	ogzstream OUTSELE ((StatSELE).c_str());
	if((!OUTSELE.good()))
	{
		cerr << "open OUT File error: "<<StatSELE<<endl;
		return 0;
	}

	double AARR;
	string chr;
	string start;
	string end;
	getline(INAA,chr);

	long double diff=0;
	long double SD=0;

	while(!INAA.eof())
	{
		string  line ;
		getline(INAA,line);
		if (line.length()<=0)  { continue ; }
		istringstream isone (line,istringstream::in);
		isone>>chr>>start>>end>>AARR;
		diff=AARR-meaWolRR;
		SD+=(diff*diff);
	}
	INAA.close();
	SD=sqrt(SD/ALL_count);


	igzstream INBB (Stat.c_str(),ifstream::in);
	if (INBB.fail())
	{
		cerr << "open Sub Group IN File error: "<<Stat<<endl;
		return  0;
	}

	getline(INBB,chr);
	double PValue, ZScore;
	OUTSS<<"#Chr\tStart\tEnd\tRLV\tSum_r^2\tCount\tZScore\tPvalue\n";
	OUTSELE<<"#Chr\tStart\tEnd\tRLV\tSum_r^2\tCount\tZScore\tPvalue\n";
	OUTSS<<"##RLV:"<<meaWolRR<<"\tSD:"<<SD<<"\tEffective windows Count:"<<ALL_count<<endl;

	while(!INBB.eof())
	{
		string  line ;
		getline(INBB,line);
		if (line.length()<=0) { continue ;}
		istringstream isone (line,istringstream::in);
		isone>>chr>>start>>end>>AARR;
		diff=AARR-meaWolRR;
		if  (diff<0)
		{
			PValue=pnorm5(AARR, meaWolRR , SD  , 1 , 0 );
		}
		else
		{
			PValue=(1-pnorm5(AARR, meaWolRR , SD  , 1 , 0 ));
		}

		ZScore=diff/SD;
		if  (ZScore<0){ZScore=0-ZScore;}
		OUTSS<<line<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::left)<<setprecision(2)<<ZScore<<"\t"<<setprecision(6)<<PValue<<endl;
		if  ( diff>0  &&  PValue<(paraFA04->Pvalue))
		{
			OUTSELE<<line<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::left)<<setprecision(2)<<ZScore<<"\t"<<setprecision(6)<<PValue<<endl;
		}
	}
	INBB.close();
	OUTSS.close();
	OUTSELE.close();

	string cc="rm  "+Stat ;
	std::system(cc.c_str()) ;

	DrawRRSliding_One(StatSS, paraFA04, para_18 );
	DrawRRDis_One( StatSS, StatSELE,  paraFA04, para_18);
	return 1;
}

#endif // processV3_H_ //

///////// swimming in the sky and flying in the sea ////////////




