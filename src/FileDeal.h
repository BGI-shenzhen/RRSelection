#ifndef ReadDataIn_H_
#define ReadDataIn_H_


using namespace std;




int Read_TwoPopVCF_IN_Phase(In3str1v * paraFA04, Para_18 * para_18 ,map <string ,int > SubVetor , map <string,map <llong, vector <BaseType> > > &  SNPListV1, map <string,map <llong, vector <BaseType> > > &  SNPListV2,  map <string,llong> & MaxSNPSite ,int & Flag_for_pro)
{

	igzstream VCFIN ((paraFA04->InStr1).c_str(),ifstream::in);
	if (VCFIN.fail())
	{
		cerr << "open VCF File IN File error: "<<(paraFA04->InStr1)<<endl;
		return  0;
	}

	vector <int> SampleSiteV1; 
	vector <int> SampleSiteV2; 
	vector<string> Vsample;
	map <string ,int >  :: iterator it;
	map <string ,int > CountSample ;
	for (it = SubVetor.begin(); it!=SubVetor.end(); it++)
	{
		CountSample.insert(map <string ,int> ::value_type(it->first,0));
	}

	while(!VCFIN.eof())
	{
		string  line ;
		getline(VCFIN,line);
		if (line.length()<=0 )  { continue ;}
		else if ( line[0] == '#' && line[1] == '#')  { continue  ; }
		else if ( line[0] == '#' && line[1] != '#')
		{
			Vsample.clear();
			split(line,Vsample," \t");
			if  ( Vsample[0]  != "#CHROM")
			{
				continue ;
			}
			int A=Vsample.size();
			for (int ii=9 ; ii< A ; ii++)
			{
				it=SubVetor.find(Vsample[ii]);
				if (it!=SubVetor.end())
				{
					if ((it->second)==1)
					{
						SampleSiteV1.push_back(ii);
					}
					else if ((it->second)==2)
					{
						SampleSiteV2.push_back(ii);
					}
					CountSample[Vsample[ii]]++;
				}
			}
			break ;
		}
		else if ( line[0] != '#' && line[1] != '#' )
		{
			cerr<<"wrong Line : "<<line<<endl;
			cerr<<"VCF Header same thing wrong, can find sample info before site info"<<endl;
			cerr<<"VCF Header sample info Flag : [  #CHROM  ] "<<endl;
			return  0;
			break;
		}
	}

	int NumberSubGroupV1=SampleSiteV1.size();
	cout<<"the Number of subPop One samples[found in VCF] is "<<NumberSubGroupV1<<endl;
	if (NumberSubGroupV1<10)
	{
		cerr<<"sub Group Population szie is too small, recommends SampleSize >=10, SubGroup One sample size: "<<NumberSubGroupV1<<endl;
		return  0;
	}

	int NumberSubGroupV2=SampleSiteV2.size();
	cout<<"the Number of subPop Tow samples[found in VCF] is "<<NumberSubGroupV2<<endl;
	if (NumberSubGroupV2<10)
	{
		cerr<<"sub Group Population szie is too small, recommends SampleSize >=10, SubGroup Two sample size: "<<NumberSubGroupV2<<endl;
		return  0;
	}


	for(it=CountSample.begin(); it!=CountSample.end() ;  it++)
	{
		if ((it->second)==0)
		{
			cerr<<"warning : Sample [ "<<(it->first)<<" ] can't be found in the VCF Header\n";
		}
		else if ((it->second)>1)
		{
			cerr<<"warning : Sample [ "<<(it->first)<<" ] can be found [Repeat] in the VCF Header\n";
		}
	}


	int Asample=Vsample.size();

	int BadSite=0;
	int BadIndelSite=0;

	vector<string> Alt ;
	vector<string> Btmp  ;
	int  Base_len ;
	llong Site ;
	int Het_count=0;
	int Miss_count=0;
	char ABase ;
	char BBase ;

	int BaseConut=0;
	char best_base='N';
	char sed_base='N';
	int Max=0;
	int SeD=0;
	map <char,int>  :: iterator itSS ;
	int ERA=0;
	vector<char>  genotype ;
	vector<BaseType>  genotypeVE ;


	while(!VCFIN.eof())
	{
		string  line ;
		getline(VCFIN,line);
		if (line.length()<=0  )  { continue  ; }
		istringstream isoneLine (line,istringstream::in);
		for (int iik=0 ; iik<Asample ; iik++)
		{
			isoneLine>>Vsample[iik];
		}

		Base_len=Vsample[3].length();
		Alt.clear();
		split(Vsample[4],Alt,",");
		for (int ii=0 ; ii<Alt.size() ;ii++)
		{
			if (Alt[ii].length()>Base_len)
			{
				Base_len=Alt[ii].length();
			}
		}

		if (Base_len>1)
		{
			BadIndelSite++;
			continue ;
		}


		for (int kkkk=0 ; kkkk<2 ;  kkkk++) 
		{
			if  (kkkk==0)
			{
				map <char,int > Count ;
				Het_count=0;
				Miss_count=0;

				for (int kk=0 ; kk< NumberSubGroupV1 ; kk++)
				{
					Btmp.clear();
					split(Vsample[SampleSiteV1[kk]], Btmp,":");
					string Genotype=Btmp[0];

					if (  Genotype[0] == '.' )
					{
						Miss_count++;
					}
					else
					{
						if (Genotype[0] != Genotype[2] )
						{
							Het_count++;
						}
						Count[Genotype[0]]++;
						Count[Genotype[2]]++;
					}
				}


				if ( ( (Miss_count*1.0/NumberSubGroupV1)  >(para_18->Miss)  )  )
				{
					continue ;
				}

				if ( ( (Het_count*1.0/NumberSubGroupV1)  >(para_18->Het) ) )
				{
					continue ;
				}

				BaseConut=0;
				best_base='N';
				sed_base='N';
				Max=0;
				SeD=0;		

				for ( itSS=Count.begin();  itSS!=Count.end(); itSS++ )
				{
					if ( (itSS->first ) == 'N' )
					{
						continue ;
					}
					else if ((itSS->second)  > Max )
					{
						SeD=Max;
						sed_base=best_base;
						Max=(itSS->second);
						best_base=itSS->first;
					}
					else if ( (itSS->second)  >= SeD )
					{
						SeD=(itSS->second);
						sed_base=itSS->first;
					}
					BaseConut++;
				}
				if (BaseConut==1 || BaseConut >2 )
				{
					BadSite++;
					continue ;
				}

				//if ( (  (1-(Max*0.5/NumberSubGroup))  < (para_18->MAF) )  )
				if ( (SeD*1.0/(SeD+Max))  < (para_18->MAF) )  
				{
					continue ;
				}

				genotypeVE.clear();
				genotype.clear();
				//		vector<BaseType>  genotypeVE ;

				BaseType  TypeA;

				//		vector<char>  genotype ;
				for (int kk=0 ; kk< NumberSubGroupV1 ; kk++)
				{
					Btmp.clear();
					split(Vsample[SampleSiteV1[kk]], Btmp,":");
					string Genotype=Btmp[0];
					ABase=Genotype[0];
					if (  ABase == '.' )
					{
						genotype.push_back('N');
						genotype.push_back('N');
					}
					else
					{
						BBase=Genotype[2];
						genotype.push_back(ABase); // phase
						genotype.push_back(BBase);
					}
				}


				ERA=genotype.size();
				for (int hh=0 ; hh<ERA ;hh++)
				{
					if (genotype[hh] == best_base )
					{
						TypeA.Value=0;
					}
					else if (genotype[hh] == sed_base )
					{
						TypeA.Value=1;
					}
					else
					{
						TypeA.Value=2; 
					}
					genotypeVE.push_back(TypeA);	
				}

				istringstream isone (Vsample[1],istringstream::in);
				isone>> Site ;


				map <string,map <llong, vector <BaseType> > >  :: iterator itSNP=SNPListV1.find(Vsample[0]);


				if (itSNP == SNPListV1.end())
				{
					map <llong, vector <BaseType> > DD;
					DD[Site]=genotypeVE;
					SNPListV1.insert(map <string,map <llong,vector <BaseType> > > ::value_type(Vsample[0],DD));
					Flag_for_pro++;
					MaxSNPSite.insert(map <string,llong>::value_type(Vsample[0],Site));
				}
				else
				{
					(itSNP->second).insert(map <llong, vector <BaseType> >  :: value_type(Site,genotypeVE));
					Flag_for_pro++;
					map <string, llong> :: iterator  MaxIT =MaxSNPSite.find(Vsample[0]);
					if ( (MaxIT->second) < Site )
					{
						MaxIT->second=Site ;
					}
				}


			}
			else 
			{





				map <char,int > Count ;
				Het_count=0;
				Miss_count=0;

				for (int kk=0 ; kk< NumberSubGroupV2 ; kk++)
				{
					Btmp.clear();
					split(Vsample[SampleSiteV2[kk]], Btmp,":");
					string Genotype=Btmp[0];

					if ( Genotype[0] == '.')
					{
						Miss_count++;
					}
					else
					{
						if (Genotype[0] != Genotype[2])
						{
							Het_count++;
						}
						Count[Genotype[0]]++;
						Count[Genotype[2]]++;
					}
				}


				if ( ( (Miss_count*1.0/NumberSubGroupV2)  >(para_18->Miss)  )  )
				{
					continue ;
				}

				if ( ( (Het_count*1.0/NumberSubGroupV2)  >(para_18->Het) ) )
				{
					continue ;
				}

				BaseConut=0;
				best_base='N';
				sed_base='N';
				Max=0;
				SeD=0;		

				for ( itSS=Count.begin();  itSS!=Count.end(); itSS++ )
				{
					if ( (itSS->first ) == 'N' )
					{
						continue ;
					}
					else if ((itSS->second)  > Max )
					{
						SeD=Max;
						sed_base=best_base;
						Max=(itSS->second);
						best_base=itSS->first;
					}
					else if ( (itSS->second)  >= SeD )
					{
						SeD=(itSS->second);
						sed_base=itSS->first;
					}
					BaseConut++;
				}
				if (BaseConut==1 || BaseConut >2 )
				{
					BadSite++;
					continue ;
				}

				//if ( (  (1-(Max*0.5/NumberSubGroup))  < (para_18->MAF) )  )
				if ( (SeD*1.0/(SeD+Max))  < (para_18->MAF) )  
				{
					continue ;
				}

				genotypeVE.clear();
				genotype.clear();
				//		vector<BaseType>  genotypeVE ;

				BaseType  TypeA;

				//		vector<char>  genotype ;
				for (int kk=0 ; kk< NumberSubGroupV2 ; kk++)
				{
					Btmp.clear();
					split(Vsample[SampleSiteV2[kk]], Btmp,":");
					string Genotype=Btmp[0];
					ABase=Genotype[0];
					if (  ABase == '.' )
					{
						genotype.push_back('N');
						genotype.push_back('N');
					}
					else
					{
						BBase=Genotype[2];
						genotype.push_back(ABase);    //  phase VCF
						genotype.push_back(BBase);    //  phase VCF
					}
				}


				ERA=genotype.size();
				for (int hh=0 ; hh<ERA ;hh++)
				{
					if (genotype[hh] == best_base )
					{
						TypeA.Value=0;
					}
					else if (genotype[hh] == sed_base )
					{
						TypeA.Value=1;
					}
					else
					{
						TypeA.Value=2; 
					}
					genotypeVE.push_back(TypeA);	
				}

				istringstream isone (Vsample[1],istringstream::in);
				isone>> Site ;


				map <string,map <llong, vector <BaseType> > >  :: iterator itSNP=SNPListV2.find(Vsample[0]);


				if (itSNP == SNPListV2.end())
				{
					map <llong, vector <BaseType> > DD;
					DD[Site]=genotypeVE;
					SNPListV2.insert(map <string,map <llong,vector <BaseType> > > ::value_type(Vsample[0],DD));
					Flag_for_pro++;
					MaxSNPSite.insert(map <string,llong>::value_type(Vsample[0],Site));
				}
				else
				{
					(itSNP->second).insert(map <llong, vector <BaseType> >  :: value_type(Site,genotypeVE));
					Flag_for_pro++;
					map <string, llong> :: iterator  MaxIT =MaxSNPSite.find(Vsample[0]);
					if ( (MaxIT->second) < Site )
					{
						MaxIT->second=Site ;
					}
				}





			}
		}

	}
	if(BadIndelSite!=0)
	{
		cout<<"warning skip Indel site, there are total skip Indel sites number is : "<<BadIndelSite<<endl;
	}

	return 1;
}







int Read_TwoPopVCF_IN(In3str1v * paraFA04, Para_18 * para_18 ,map <string ,int > SubVetor , map <string,map <llong, vector <BaseType> > > &  SNPListV1, map <string,map <llong, vector <BaseType> > > &  SNPListV2,  map <string,llong> & MaxSNPSite ,int & Flag_for_pro)
{

	igzstream VCFIN ((paraFA04->InStr1).c_str(),ifstream::in);
	if (VCFIN.fail())
	{
		cerr << "open VCF File IN File error: "<<(paraFA04->InStr1)<<endl;
		return  0;
	}

	vector <int> SampleSiteV1; 
	vector <int> SampleSiteV2; 
	vector<string> Vsample ;
	map <string ,int >  :: iterator it;
	map <string ,int > CountSample;
	for (it = SubVetor.begin(); it!=SubVetor.end(); it++)
	{
		CountSample.insert(map <string ,int> ::value_type(it->first,0));
	}

	while(!VCFIN.eof())
	{
		string  line ;
		getline(VCFIN,line);
		if (line.length()<=0 )  { continue ;}
		else if ( line[0] == '#' && line[1] == '#')  { continue  ; }
		else if ( line[0] == '#' && line[1] != '#')
		{
			Vsample.clear();
			split(line,Vsample," \t");
			if  ( Vsample[0]  != "#CHROM")
			{
				continue ;
			}
			int A=Vsample.size();
			for (int ii=9 ; ii< A ; ii++)
			{
				it=SubVetor.find(Vsample[ii]);
				if (it!=SubVetor.end())
				{
					if ((it->second)==1)
					{
						SampleSiteV1.push_back(ii);
					}
					else if ((it->second)==2)
					{
						SampleSiteV2.push_back(ii);
					}
					CountSample[Vsample[ii]]++;
				}
			}
			break ;
		}
		else if ( line[0] != '#' && line[1] != '#' )
		{
			cerr<<"wrong Line : "<<line<<endl;
			cerr<<"VCF Header same thing wrong, can find sample info before site info"<<endl;
			cerr<<"VCF Header sample info Flag : [  #CHROM  ] "<<endl;
			return  0;
			break;
		}
	}

	int NumberSubGroupV1=SampleSiteV1.size();
	cout<<"the Number of subPop One samples[found in VCF] is "<<NumberSubGroupV1<<endl;
	if (NumberSubGroupV1<10)
	{
		cerr<<"sub Group Population szie is too small, recommends SampleSize >=10, SubGroup One sample size: "<<NumberSubGroupV1<<endl;
		return  0;
	}

	int NumberSubGroupV2=SampleSiteV2.size();
	cout<<"the Number of subPop Tow samples[found in VCF] is "<<NumberSubGroupV2<<endl;
	if (NumberSubGroupV2<10)
	{
		cerr<<"sub Group Population szie is too small, recommends SampleSize >=10, SubGroup Two sample size: "<<NumberSubGroupV2<<endl;
		return  0;
	}


	for(it=CountSample.begin(); it!=CountSample.end() ;  it++)
	{
		if ((it->second)==0)
		{
			cerr<<"warning : Sample [ "<<(it->first)<<" ] can't be found in the VCF Header\n";
		}
		else if ((it->second)>1)
		{
			cerr<<"warning : Sample [ "<<(it->first)<<" ] can be found [Repeat] in the VCF Header\n";
		}
	}


	int Asample=Vsample.size();

	int BadSite=0;
	int BadIndelSite=0;

	vector<string> Alt ;
	vector<string> Btmp  ;
	int  Base_len ;
	llong Site ;
	int Het_count=0;
	int Miss_count=0;
	char ABase ;
	char BBase ;

	int BaseConut=0;
	char best_base='N';
	char sed_base='N';
	int Max=0;
	int SeD=0;
	map <char,int>  :: iterator itSS ;
	int ERA=0;
	vector<char>  genotype ;
	vector<BaseType>  genotypeVE ;


	while(!VCFIN.eof())
	{
		string  line ;
		getline(VCFIN,line);
		if (line.length()<=0  )  { continue  ; }
		istringstream isoneLine (line,istringstream::in);
		for (int iik=0 ; iik<Asample ; iik++)
		{
			isoneLine>>Vsample[iik];
		}

		Base_len=Vsample[3].length();
		Alt.clear();
		split(Vsample[4],Alt,",");
		for (int ii=0 ; ii<Alt.size() ;ii++)
		{
			if (Alt[ii].length()>Base_len)
			{
				Base_len=Alt[ii].length();
			}
		}

		if (Base_len>1)
		{
			BadIndelSite++;
			continue ;
		}


		for (int kkkk=0 ; kkkk<2 ;  kkkk++) 
		{
			if  (kkkk==0)
			{
				map <char,int > Count ;
				Het_count=0;
				Miss_count=0;

				for (int kk=0 ; kk< NumberSubGroupV1 ; kk++)
				{
					Btmp.clear();
					split(Vsample[SampleSiteV1[kk]], Btmp,":");
					string Genotype=Btmp[0];

					if (  Genotype[0] == '.' )
					{
						Miss_count++;
					}
					else
					{
						if (Genotype[0] != Genotype[2] )
						{
							Het_count++;
						}
						Count[Genotype[0]]++;
						Count[Genotype[2]]++;
					}
				}


				if ( ( (Miss_count*1.0/NumberSubGroupV1)  >(para_18->Miss)  )  )
				{
					continue ;
				}

				if ( ( (Het_count*1.0/NumberSubGroupV1)  >(para_18->Het) ) )
				{
					continue ;
				}

				BaseConut=0;
				best_base='N';
				sed_base='N';
				Max=0;
				SeD=0;		

				for ( itSS=Count.begin();  itSS!=Count.end(); itSS++ )
				{
					if ( (itSS->first ) == 'N' )
					{
						continue ;
					}
					else if ((itSS->second)  > Max )
					{
						SeD=Max;
						sed_base=best_base;
						Max=(itSS->second);
						best_base=itSS->first;
					}
					else if ( (itSS->second)  >= SeD )
					{
						SeD=(itSS->second);
						sed_base=itSS->first;
					}
					BaseConut++;
				}
				if (BaseConut==1 || BaseConut >2 )
				{
					BadSite++;
					continue ;
				}

				//if ( (  (1-(Max*0.5/NumberSubGroup))  < (para_18->MAF) )  )
				if ( (SeD*1.0/(SeD+Max))  < (para_18->MAF) )  
				{
					continue ;
				}

				genotypeVE.clear();
				genotype.clear();
				//		vector<BaseType>  genotypeVE ;

				BaseType  TypeA;

				//		vector<char>  genotype ;
				for (int kk=0 ; kk< NumberSubGroupV1 ; kk++)
				{
					Btmp.clear();
					split(Vsample[SampleSiteV1[kk]], Btmp,":");
					string Genotype=Btmp[0];
					ABase=Genotype[0];
					if (  ABase == '.' )
					{
						genotype.push_back('N');
						genotype.push_back('N');
					}
					else
					{
						BBase=Genotype[2];
						if  (ABase != BBase)
						{
							genotype.push_back(best_base);   // best and sed base to phase
							genotype.push_back(sed_base);
						}
						else
						{
							genotype.push_back(ABase);
							genotype.push_back(BBase);
						}
					}
				}


				ERA=genotype.size();
				for (int hh=0 ; hh<ERA ;hh++)
				{
					if (genotype[hh] == best_base )
					{
						TypeA.Value=0;
					}
					else if (genotype[hh] == sed_base )
					{
						TypeA.Value=1;
					}
					else
					{
						TypeA.Value=2; 
					}
					genotypeVE.push_back(TypeA);	
				}

				istringstream isone (Vsample[1],istringstream::in);
				isone>> Site ;


				map <string,map <llong, vector <BaseType> > >  :: iterator itSNP=SNPListV1.find(Vsample[0]);


				if (itSNP == SNPListV1.end())
				{
					map <llong, vector <BaseType> > DD;
					DD[Site]=genotypeVE;
					SNPListV1.insert(map <string,map <llong,vector <BaseType> > > ::value_type(Vsample[0],DD));
					Flag_for_pro++;
					MaxSNPSite.insert(map <string,llong>::value_type(Vsample[0],Site));
				}
				else
				{
					(itSNP->second).insert(map <llong, vector <BaseType> >  :: value_type(Site,genotypeVE));
					Flag_for_pro++;
					map <string, llong> :: iterator  MaxIT =MaxSNPSite.find(Vsample[0]);
					if ( (MaxIT->second) < Site )
					{
						MaxIT->second=Site ;
					}
				}


			}
			else 
			{





				map <char,int > Count ;
				Het_count=0;
				Miss_count=0;

				for (int kk=0 ; kk< NumberSubGroupV2 ; kk++)
				{
					Btmp.clear();
					split(Vsample[SampleSiteV2[kk]], Btmp,":");
					string Genotype=Btmp[0];

					if ( Genotype[0] == '.')
					{
						Miss_count++;
					}
					else
					{
						if (Genotype[0] != Genotype[2])
						{
							Het_count++;
						}
						Count[Genotype[0]]++;
						Count[Genotype[2]]++;
					}
				}


				if ( ( (Miss_count*1.0/NumberSubGroupV2)  >(para_18->Miss)  )  )
				{
					continue ;
				}

				if ( ( (Het_count*1.0/NumberSubGroupV2)  >(para_18->Het) ) )
				{
					continue ;
				}

				BaseConut=0;
				best_base='N';
				sed_base='N';
				Max=0;
				SeD=0;		

				for ( itSS=Count.begin();  itSS!=Count.end(); itSS++ )
				{
					if ( (itSS->first ) == 'N' )
					{
						continue ;
					}
					else if ((itSS->second)  > Max )
					{
						SeD=Max;
						sed_base=best_base;
						Max=(itSS->second);
						best_base=itSS->first;
					}
					else if ( (itSS->second)  >= SeD )
					{
						SeD=(itSS->second);
						sed_base=itSS->first;
					}
					BaseConut++;
				}
				if (BaseConut==1 || BaseConut >2 )
				{
					BadSite++;
					continue ;
				}

				//if ( (  (1-(Max*0.5/NumberSubGroup))  < (para_18->MAF) )  )
				if ( (SeD*1.0/(SeD+Max))  < (para_18->MAF) )  
				{
					continue ;
				}

				genotypeVE.clear();
				genotype.clear();
				//		vector<BaseType>  genotypeVE ;

				BaseType  TypeA;

				//		vector<char>  genotype ;
				for (int kk=0 ; kk< NumberSubGroupV2 ; kk++)
				{
					Btmp.clear();
					split(Vsample[SampleSiteV2[kk]], Btmp,":");
					string Genotype=Btmp[0];
					ABase=Genotype[0];
					if (  ABase == '.' )
					{
						genotype.push_back('N');
						genotype.push_back('N');
					}
					else
					{
						BBase=Genotype[2];
						if  (ABase != BBase)
						{
							genotype.push_back(best_base);   // best and sed base to phase
							genotype.push_back(sed_base);
						}
						else
						{
							genotype.push_back(ABase);
							genotype.push_back(BBase);
						}
					}
				}


				ERA=genotype.size();
				for (int hh=0 ; hh<ERA ;hh++)
				{
					if (genotype[hh] == best_base )
					{
						TypeA.Value=0;
					}
					else if (genotype[hh] == sed_base )
					{
						TypeA.Value=1;
					}
					else
					{
						TypeA.Value=2; 
					}
					genotypeVE.push_back(TypeA);	
				}

				istringstream isone (Vsample[1],istringstream::in);
				isone>> Site ;


				map <string,map <llong, vector <BaseType> > >  :: iterator itSNP=SNPListV2.find(Vsample[0]);


				if (itSNP == SNPListV2.end())
				{
					map <llong, vector <BaseType> > DD;
					DD[Site]=genotypeVE;
					SNPListV2.insert(map <string,map <llong,vector <BaseType> > > ::value_type(Vsample[0],DD));
					Flag_for_pro++;
					MaxSNPSite.insert(map <string,llong>::value_type(Vsample[0],Site));
				}
				else
				{
					(itSNP->second).insert(map <llong, vector <BaseType> >  :: value_type(Site,genotypeVE));
					Flag_for_pro++;
					map <string, llong> :: iterator  MaxIT =MaxSNPSite.find(Vsample[0]);
					if ( (MaxIT->second) < Site )
					{
						MaxIT->second=Site ;
					}
				}





			}
		}

	}
	if(BadIndelSite!=0)
	{
		cout<<"warning skip Indel site, there are total skip Indel sites number is : "<<BadIndelSite<<endl;
	}

	return 1;
}




int Read_SubPopVCF_IN(In3str1v * paraFA04, Para_18 * para_18 , map <string,map <llong, vector <BaseType> > > &  SNPList, map <string,llong> & MaxSNPSite ,int & Flag_for_pro)
{
	igzstream SampleList ((paraFA04->SubPop).c_str(),ifstream::in);
	if (SampleList.fail())
	{
		cerr << "open Sub Group IN File error: "<<(paraFA04->SubPop)<<endl;
		return 0;
	}

	map <string ,int >  SubVetor;
	map <string ,int >  :: iterator it;

	while(!SampleList.eof())
	{
		string  line ;
		getline(SampleList,line);
		if (line.length()<=0 || line[0] == '#' )  { continue  ; }
		line =line.substr(line.find(':') ==string::npos ? 0 : line.find(':') + 1);
		vector<string> inf;
		split(line,inf," \t");
		int A=inf.size();
		for(int ii=0 ; ii<A; ii++)
		{
			it=SubVetor.find(inf[ii]);
			if (it==SubVetor.end())
			{
				SubVetor.insert(map <string ,int> ::value_type(inf[ii],0));
			}
		}
	}
	SampleList.close();

	igzstream VCFIN ((paraFA04->InStr1).c_str(),ifstream::in);
	if (VCFIN.fail())
	{
		cerr << "open VCF File IN File error: "<<(paraFA04->InStr1)<<endl;
		return  0;
	}

	vector <int> SampleSite; 
	vector<string> Vsample ;

	while(!VCFIN.eof())
	{
		string  line ;
		getline(VCFIN,line);
		if (line.length()<=0 )  { continue  ; }
		else if ( line[0] == '#' && line[1] == '#' )  { continue  ; }
		else if ( line[0] == '#' && line[1] != '#' )
		{
			Vsample.clear();
			split(line,Vsample," \t");
			if  ( Vsample[0]  != "#CHROM")
			{
				continue  ;
			}
			int A=Vsample.size();

			for (int ii=9 ; ii< A ; ii++)
			{
				it=SubVetor.find(Vsample[ii]);
				if (it!=SubVetor.end())
				{
					SampleSite.push_back(ii);
					(it->second)++;
				}
			}
			break ;
		}
		else if ( line[0] != '#' && line[1] != '#' )
		{
			cerr<<"wrong Line : "<<line<<endl;
			cerr<<"VCF Header same thing wrong, can find sample info before site info"<<endl;
			cerr<<"VCF Header sample info Flag : [  #CHROM  ] "<<endl;
			return  0;
			break;
		}
	}

	int NumberSubGroup=SampleSite.size();
	cout<<"the Number of subPop samples[found in VCF] is "<<NumberSubGroup<<endl;
	if (NumberSubGroup<10)
	{
		cerr<<"sub Group Population szie is too small, recommends SampleSize >=10, SubGroup sample size: "<<NumberSubGroup<<endl;
		return  0;
	}


	for(it=SubVetor.begin(); it!=SubVetor.end() ;  it++)
	{
		if ((it->second)==0)
		{
			cerr<<"warning : Sample [ "<<(it->first)<<" ] can't be found in the VCF Header\n";
		}
		else if ((it->second)>1)
		{
			cerr<<"warning : Sample [ "<<(it->first)<<" ] can be found [Repeat] in the VCF Header\n";
		}
	}

	int Asample=Vsample.size();

	int BadSite=0;
	int BadIndelSite=0;

	vector<string> Alt ;
	vector<string> Btmp  ;
	int  Base_len ;
	llong Site ;
	int Het_count=0;
	int Miss_count=0;
	char ABase ;
	char BBase ;

	int BaseConut=0;
	char best_base='N';
	char sed_base='N';
	int Max=0;
	int SeD=0;
	map <char,int>  :: iterator itSS ;
	int ERA=0;
	vector<char>  genotype ;
	vector<BaseType>  genotypeVE ;


	while(!VCFIN.eof())
	{
		string  line ;
		getline(VCFIN,line);
		if (line.length()<=0  )  { continue  ; }
		istringstream isoneLine (line,istringstream::in);
		for (int iik=0 ; iik<Asample ; iik++)
		{
			isoneLine>>Vsample[iik];
		}

		Base_len=Vsample[3].length();
		Alt.clear();
		split(Vsample[4],Alt,",");
		for (int ii=0 ; ii<Alt.size() ;ii++)
		{
			if (Alt[ii].length()>Base_len)
			{
				Base_len=Alt[ii].length();
			}
		}

		if (Base_len>1)
		{
			BadIndelSite++;
			continue ;
		}

		map <char,int > Count ;
		Het_count=0;
		Miss_count=0;

		for (int kk=0 ; kk< NumberSubGroup ; kk++)
		{
			Btmp.clear();
			split(Vsample[SampleSite[kk]], Btmp,":");
			string Genotype=Btmp[0];

			if (  Genotype[0] == '.' )
			{
				Miss_count++ ;
			}
			else
			{
				if  (Genotype[0] != Genotype[2] )
				{
					Het_count++;
				}
				Count[Genotype[0]]++;
				Count[Genotype[2]]++;
			}
		}


		if ( ( (Miss_count*1.0/NumberSubGroup)  >(para_18->Miss)  )  )
		{
			continue ;
		}

		if ( ( (Het_count*1.0/NumberSubGroup)  >(para_18->Het) )  )
		{
			continue ;
		}

		BaseConut=0;
		best_base='N';
		sed_base='N';
		Max=0;
		SeD=0;		

		for ( itSS=Count.begin();  itSS!=Count.end(); itSS++ )
		{
			if ( (itSS->first ) == 'N' )
			{
				continue ;
			}
			else if ((itSS->second)  > Max )
			{
				SeD=Max;
				sed_base=best_base;
				Max=(itSS->second);
				best_base=itSS->first;
			}
			else if ( (itSS->second)  >= SeD )
			{
				SeD=(itSS->second);
				sed_base=itSS->first;
			}
			BaseConut++;
		}
		if (BaseConut==1 || BaseConut >2 )
		{
			BadSite++;
			continue ;
		}

		//if ( (  (1-(Max*0.5/NumberSubGroup))  < (para_18->MAF) )  )
		if ( (SeD*1.0/(SeD+Max))  < (para_18->MAF) )  
		{
			continue ;
		}

		genotypeVE.clear();
		genotype.clear();
		//		vector<BaseType>  genotypeVE ;

		BaseType  TypeA;

		//		vector<char>  genotype ;
		for (int kk=0 ; kk< NumberSubGroup ; kk++)
		{
			Btmp.clear();
			split(Vsample[SampleSite[kk]], Btmp,":");
			string Genotype=Btmp[0];
			ABase=Genotype[0];
			if (  ABase == '.' )
			{
				genotype.push_back('N');
				genotype.push_back('N');
			}
			else
			{
				BBase=Genotype[2];
				if  (ABase != BBase)
				{
					genotype.push_back(best_base);   // best and sed base to phase
					genotype.push_back(sed_base);
				}
				else
				{
					genotype.push_back(ABase);
					genotype.push_back(BBase);
				}
			}
		}


		ERA=genotype.size();
		for (int hh=0 ; hh<ERA ;hh++)
		{
			if (genotype[hh] == best_base )
			{
				TypeA.Value=0;
			}
			else if (genotype[hh] == sed_base )
			{
				TypeA.Value=1;
			}
			else
			{
				TypeA.Value=2; 
			}
			genotypeVE.push_back(TypeA);	
		}

		istringstream isone (Vsample[1],istringstream::in);
		isone>> Site ;


		map <string,map <llong, vector <BaseType> > >  :: iterator itSNP=SNPList.find(Vsample[0]);


		if (itSNP == SNPList.end())
		{
			map <llong, vector <BaseType> > DD;
			DD[Site]=genotypeVE;
			SNPList.insert(map <string,map <llong,vector <BaseType> > > ::value_type(Vsample[0],DD));
			Flag_for_pro++;
			MaxSNPSite.insert(map <string,llong>::value_type(Vsample[0],Site));
		}
		else
		{
			(itSNP->second).insert(map <llong, vector <BaseType> >  :: value_type(Site,genotypeVE));
			Flag_for_pro++;
			map <string, llong> :: iterator  MaxIT =MaxSNPSite.find(Vsample[0]);
			if ( (MaxIT->second) < Site )
			{
				MaxIT->second=Site ;
			}

		}
	}

	if(BadIndelSite!=0)
	{
		cout<<"warning skip Indel site, there are total skip Indel sites number is : "<<BadIndelSite<<endl;
	}
	if (BadSite!=0)
	{
		cout<<"Warning skip non bi-allelic(Singleton/ThreeMulti allelic) site, and total skip allelic sites number is :"<<BadSite<<endl;
	}

	return 1;
}



bool DeleVCFPhase( string  VCFfile )
{

	igzstream VCFIN ( VCFfile.c_str(),ifstream::in);
	if (VCFIN.fail())
	{
		cerr << "open VCF File IN File error: "<<VCFfile<<endl;
		return  false;
	}

	vector<string> inf ;
	bool  TTFF=false ;
	while(!VCFIN.eof())
	{
		string  line ;
		getline(VCFIN,line);
		if (line.length()<=0 )  { continue  ; }
		else if ( line[0] == '#'  )  { continue  ; }
		else
		{
			split(line,inf," \t");
			if ((inf[9])[1] == '|') 
			{
				TTFF=true;
			}
			break ;
		}
	}
	VCFIN.close();
	if  (TTFF)
	{
		cout <<"#Detected VCF File is phased file with '|', Read VCF in Phase mode"<<endl;
	}
	return TTFF;
}






int Read_VCF_IN(In3str1v * paraFA04, Para_18 * para_18 , map <string,map <llong, vector <BaseType> > > &  SNPList ,map <string,llong> & MaxSNPSite,int & Flag_for_pro )
{

	igzstream VCFIN ((paraFA04->InStr1).c_str(),ifstream::in);

	if (VCFIN.fail())
	{
		cerr << "open VCF File IN File error: "<<(paraFA04->InStr1)<<endl;
		return  0;
	}

	vector<string> inf ;

	while(!VCFIN.eof())
	{
		string  line ;
		getline(VCFIN,line);
		if (line.length()<=0 )  { continue  ; }
		else if ( line[0] == '#' && line[1] == '#' )  { continue  ; }
		else if ( line[0] == '#' && line[1] != '#' )
		{
			inf.clear();
			split(line,inf," \t");
			if  ( inf[0]  != "#CHROM")
			{
				continue  ;
			}			
			break ;
		}
		else if ( line[0] != '#' && line[1] != '#' )
		{
			cerr<<"wrong Line : "<<line<<endl;
			cerr<<"VCF Header same thing wrong, can find sample info before site info"<<endl;
			cerr<<"VCF Header sample info Flag : [  #CHROM  ] "<<endl;
			return  0;
			break;
		}
	}

	int BadSite=0 ;
	int BadIndelSite=0;

	int Asample=inf.size();
	int SampleNum=(Asample-9);

	if (SampleNum<10)
	{
		cerr<<"sub Group Population szie is too small, recommends SampleSize >=10,but sample size: "<<SampleNum<<endl;
		return  0;
	}



	vector<string> Alt ;
	vector<string> Btmp  ;
	int  Base_len ;
	llong Site ;
	int Het_count=0;
	int Miss_count=0;

	int BaseConut=0;
	char best_base='N';
	char sed_base='N';
	int Max=0;
	int SeD=0;
	map <char,int>  :: iterator it ;
	int ERA=0;
	vector<char>  genotype ;
	vector<BaseType>  genotypeVE ;



	while(!VCFIN.eof())
	{
		string  line ;
		getline(VCFIN,line);
		if (line.length()<=0)  { continue  ; }

		istringstream isoneLine (line,istringstream::in);
		for (int iik=0 ; iik<Asample ; iik++)
		{
			isoneLine>>inf[iik];
		}
		Base_len=inf[3].length();
		Alt.clear();
		split(inf[4],Alt,",");

		for (int ii=0 ; ii<Alt.size() ;ii++)
		{
			if (Alt[ii].length()>Base_len)
			{
				Base_len=Alt[ii].length();
			}
		}

		if (Base_len>1)
		{
			BadIndelSite++;
			continue ;
		}

		map <char,int > Count ;
		Het_count=0;
		Miss_count=0;

		for (int jj=9 ; jj< Asample ;jj++ )
		{
			Btmp.clear();
			split(inf[jj], Btmp,":");
			string Genotype=Btmp[0];
			if (  Genotype[0] == '.' )
			{
				Miss_count++ ;
			}
			else
			{
				if  (Genotype[0] != Genotype[2] )
				{
					Het_count++;
				}
				Count[Genotype[0]]++;
				Count[Genotype[2]]++;
			}
		}

		if ( ( (Miss_count*1.0/SampleNum)  >(para_18->Miss)  )  )
		{
			continue ;
		}

		if ( ( (Het_count*1.0/SampleNum)  >(para_18->Het) )  )
		{
			continue ;
		}

		BaseConut=0;
		best_base='N';
		sed_base='N';
		Max=0;
		SeD=0;		

		for ( it=Count.begin(); it!=Count.end(); it++ )
		{
			if ( (it->first ) == 'N' )
			{					
				continue ;
			}
			else if ((it->second)  > Max )
			{
				SeD=Max;
				sed_base=best_base;
				Max=(it->second);
				best_base=it->first;
			}
			else if ( (it->second)  >= SeD )
			{
				SeD=(it->second);
				sed_base=it->first;
			}
			BaseConut++;
		}
		if (BaseConut==1 || BaseConut >2 )
		{
			BadSite++;
			continue ;
		}

		if ( (SeD*1.0/(SeD+Max))  < (para_18->MAF) )  
		{
			continue ;
		}

		genotype.clear();

		for (int jj=9 ; jj< Asample ;jj++ )
		{
			Btmp.clear();
			split(inf[jj], Btmp,":");
			string Genotype=Btmp[0];
			if (  Genotype[0] == '.' )
			{
				genotype.push_back('N');
				genotype.push_back('N');
			}
			else
			{
				if  (Genotype[0] != Genotype[2] )	
				{
					genotype.push_back(best_base);
					genotype.push_back(sed_base);
				}	
				else
				{
					genotype.push_back(Genotype[0]);
					genotype.push_back(Genotype[2]);
				}
			}
		}



		genotypeVE.clear();

		BaseType  TypeA;

		ERA=genotype.size();

		for (int hh=0 ; hh<ERA ;hh++)
		{

			if (genotype[hh] == best_base  )
			{
				TypeA.Value=0;
			}
			else if (genotype[hh] == sed_base )
			{
				TypeA.Value=1;
			}
			else
			{
				TypeA.Value=2; 
			}
			genotypeVE.push_back(TypeA);
		}





		istringstream isone (inf[1],istringstream::in);
		isone>> Site ;


		map <string,map <llong, vector <BaseType> > >  :: iterator itSNP=SNPList.find(inf[0]);

		if (itSNP == SNPList.end())
		{
			map <llong, vector <BaseType> > DD;
			DD[Site]=genotypeVE;
			SNPList.insert(map <string,map <llong,vector <BaseType> > > ::value_type(inf[0],DD));
			Flag_for_pro++;
			MaxSNPSite.insert(map <string,llong>::value_type(inf[0],Site));
		}
		else
		{
			(itSNP->second).insert(map <llong, vector <BaseType> >  :: value_type(Site,genotypeVE)) ;
			Flag_for_pro++;
			map <string, llong> :: iterator  MaxIT =MaxSNPSite.find(inf[0]);
			if ( (MaxIT->second) < Site )
			{
				MaxIT->second=Site ;
			}

		}
	}



	VCFIN.close();

	if(BadIndelSite!=0)
	{
		cout<<"warning skip Indel site, there are total skip Indel sites number is : "<<BadIndelSite<<endl;
	}
	if (BadSite!=0)
	{
		cout<<"Warning skip non bi-allelic(Singleton/ThreeMulti allelic) site, and total skip allelic sites number is :"<<BadSite<<endl;
	}

	return 1;
}





/*
inline int OUTStatFile( In3str1v * paraFA04, Para_18 * para_18 , StarRsult *All_Stat )
{
	string Stat=(paraFA04->InStr2)+".stat.gz";
	ogzstream OUT ((Stat).c_str());	
	if((!OUT.good()))
	{
		cerr << "open OUT File error: "<<Stat<<endl;
		return  0;
	}

	OUT<<"#Dist\tMean_r^2\tMean_D'\tSum_r^2\tSum_D'\tNumberPairs\n";
	(paraFA04->WindowSize)++;

	for (int ii=1 ; ii<(paraFA04->WindowSize) ; ii++ )
	{
		int count=(All_Stat[ii]).Count;
		if  (count==0)	{continue ;}
		double SumRR=(All_Stat[ii]).sumRR;
		double SumD=(All_Stat[ii]).sumD;

		double MeanRR=SumRR/count;
		double MeanD=SumD/count;
		OUT<<ii<<setprecision(4)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<"\t"<<MeanRR<<"\t"<<MeanD<<"\t"<<SumRR<<"\t"<<SumD<<"\t"<<count<<"\n";
	}

	OUT.close();
	return 1  ;
}
*/





/////////////////////   Phase  VCF /////////////////////


int Read_SubPopVCF_IN_Phase(In3str1v * paraFA04, Para_18 * para_18 , map <string,map <llong, vector <BaseType> > > &  SNPList, map <string,llong> & MaxSNPSite,int & Flag_for_pro )
{
	igzstream SampleList ((paraFA04->SubPop).c_str(),ifstream::in);
	if (SampleList.fail())
	{
		cerr << "open Sub Group IN File error: "<<(paraFA04->SubPop)<<endl;
		return  0;
	}

	map <string ,int >  SubVetor;
	map <string ,int >  :: iterator it;

	while(!SampleList.eof())
	{
		string  line ;
		getline(SampleList,line);
		if (line.length()<=0 || line[0] == '#' )  { continue  ; }
		line =line.substr(line.find(':') ==string::npos ? 0 : line.find(':') + 1);
		vector<string> inf;
		split(line,inf," \t");
		int A=inf.size();
		for(int ii=0 ; ii<A; ii++)
		{
			it=SubVetor.find(inf[ii]);
			if (it==SubVetor.end())
			{
				SubVetor.insert(map <string ,int> ::value_type(inf[ii],0));
			}
		}
	}
	SampleList.close();

	igzstream VCFIN ((paraFA04->InStr1).c_str(),ifstream::in);
	if (VCFIN.fail())
	{
		cerr << "open VCF File IN File error: "<<(paraFA04->InStr1)<<endl;
		return  0;
	}

	vector <int> SampleSite; 
	vector<string> Vsample ;

	while(!VCFIN.eof())
	{
		string  line ;
		getline(VCFIN,line);
		if (line.length()<=0 )  { continue  ; }
		else if ( line[0] == '#' && line[1] == '#' )  { continue  ; }
		else if ( line[0] == '#' && line[1] != '#' )
		{
			Vsample.clear();
			split(line,Vsample," \t");
			if  ( Vsample[0]  != "#CHROM")
			{
				continue  ;
			}
			int A=Vsample.size();

			for (int ii=9 ; ii< A ; ii++)
			{
				it=SubVetor.find(Vsample[ii]);
				if (it!=SubVetor.end())
				{
					SampleSite.push_back(ii);
					(it->second)++;
				}
			}
			break ;
		}
		else if ( line[0] != '#' && line[1] != '#' )
		{
			cerr<<"wrong Line : "<<line<<endl;
			cerr<<"VCF Header same thing wrong, can find sample info before site info"<<endl;
			cerr<<"VCF Header sample info Flag : [  #CHROM  ] "<<endl;
			return  0;
			break;
		}
	}

	int NumberSubGroup=SampleSite.size();
	cout<<"the Number of subPop samples[found in VCF] is "<<NumberSubGroup<<endl;
	if (NumberSubGroup<10)
	{
		cerr<<"sub Group Population szie is too small, recommends SampleSize >=10, SubGroup sample size: "<<NumberSubGroup<<endl;
		return  0;
	}


	for(it=SubVetor.begin(); it!=SubVetor.end() ;  it++)
	{
		if ((it->second)==0)
		{
			cerr<<"warning : Sample [ "<<(it->first)<<" ] can't be found in the VCF Header\n";
		}
		else if ((it->second)>1)
		{
			cerr<<"warning : Sample [ "<<(it->first)<<" ] can be found [Repeat] in the VCF Header\n";
		}
	}


	int Asample=Vsample.size();

	int BadSite=0;
	int BadIndelSite=0;

	vector<string> Alt ;
	vector<string> Btmp  ;
	int  Base_len ;
	llong Site ;
	int Het_count=0;
	int Miss_count=0;
	char ABase ;
	char BBase ;

	int BaseConut=0;
	char best_base='N';
	char sed_base='N';
	int Max=0;
	int SeD=0;
	map <char,int>  :: iterator itSS ;
	int ERA=0;
	vector<char>  genotype ;
	vector<BaseType>  genotypeVE ;


	while(!VCFIN.eof())
	{
		string  line ;
		getline(VCFIN,line);
		if (line.length()<=0  )  { continue  ; }
		istringstream isoneLine (line,istringstream::in);
		for (int iik=0 ; iik<Asample ; iik++)
		{
			isoneLine>>Vsample[iik];
		}

		Base_len=Vsample[3].length();
		Alt.clear();
		split(Vsample[4],Alt,",");
		for (int ii=0 ; ii<Alt.size() ;ii++)
		{
			if (Alt[ii].length()>Base_len)
			{
				Base_len=Alt[ii].length();
			}
		}

		if (Base_len>1)
		{
			BadIndelSite++;
			continue ;
		}

		map <char,int > Count ;
		Het_count=0;
		Miss_count=0;

		for (int kk=0 ; kk< NumberSubGroup ; kk++)
		{
			Btmp.clear();
			split(Vsample[SampleSite[kk]], Btmp,":");
			string Genotype=Btmp[0];

			if (  Genotype[0] == '.' )
			{
				Miss_count++ ;
			}
			else
			{
				if  (Genotype[0] != Genotype[2] )
				{
					Het_count++;
				}
				Count[Genotype[0]]++;
				Count[Genotype[2]]++;
			}
		}


		if ( ( (Miss_count*1.0/NumberSubGroup)  >(para_18->Miss)  )  )
		{
			continue ;
		}

		if ( ( (Het_count*1.0/NumberSubGroup)  >(para_18->Het) )  )
		{
			continue ;
		}

		BaseConut=0;
		best_base='N';
		sed_base='N';
		Max=0;
		SeD=0;		

		for ( itSS=Count.begin();  itSS!=Count.end(); itSS++ )
		{
			if ( (itSS->first ) == 'N' )
			{
				continue ;
			}
			else if ((itSS->second)  > Max )
			{
				SeD=Max;
				sed_base=best_base;
				Max=(itSS->second);
				best_base=itSS->first;
			}
			else if ( (itSS->second)  >= SeD )
			{
				SeD=(itSS->second);
				sed_base=itSS->first;
			}
			BaseConut++;
		}

		if (BaseConut==1 || BaseConut >2 )
		{
			BadSite++;
			continue ;
		}

		//if ( (  (1-(Max*0.5/NumberSubGroup))  < (para_18->MAF) )  )
		if ( (SeD*1.0/(SeD+Max))  < (para_18->MAF) )  
		{
			continue ;
		}

		genotypeVE.clear();
		genotype.clear();


		BaseType  TypeA;

		for (int kk=0 ; kk< NumberSubGroup ; kk++)
		{
			Btmp.clear();
			split(Vsample[SampleSite[kk]], Btmp,":");
			string Genotype=Btmp[0];
			ABase=Genotype[0];
			if (  ABase == '.' )
			{
				genotype.push_back('N');
				genotype.push_back('N');
			}
			else
			{
				BBase=Genotype[2];
				genotype.push_back(ABase);	// phase VCF
				genotype.push_back(BBase);  // phase VCF
			}
		}


		ERA=genotype.size();
		for (int hh=0 ; hh<ERA ;hh++)
		{
			if (genotype[hh] == best_base )
			{
				TypeA.Value=0;
			}
			else if (genotype[hh] == sed_base )
			{
				TypeA.Value=1;
			}
			else
			{
				TypeA.Value=2;
			}
			genotypeVE.push_back(TypeA);
		}



		istringstream isone (Vsample[1],istringstream::in);
		isone>> Site ;


		map <string,map <llong, vector <BaseType> > >  :: iterator itSNP=SNPList.find(Vsample[0]);

		if (itSNP == SNPList.end())
		{
			map <llong, vector <BaseType> > DD;
			DD[Site]=genotypeVE;
			SNPList.insert(map <string,map <llong,vector <BaseType> > > ::value_type(Vsample[0],DD));
			Flag_for_pro++;
			MaxSNPSite.insert(map <string,llong>::value_type(Vsample[0],Site));
		}
		else
		{
			(itSNP->second).insert(map <llong, vector <BaseType> >  :: value_type(Site,genotypeVE)) ;
			Flag_for_pro++;
			map <string, llong> :: iterator  MaxIT =MaxSNPSite.find(Vsample[0]);
			if ( (MaxIT->second) < Site )
			{
				MaxIT->second=Site ;
			}

		}
	}

	if(BadIndelSite!=0)
	{
		cout<<"warning skip Indel site, there are total skip Indel sites number is : "<<BadIndelSite<<endl;
	}
	if (BadSite!=0)
	{
		cout<<"Warning skip non bi-allelic(Singleton/ThreeMulti allelic) site, and total skip allelic sites number is :"<<BadSite<<endl;
	}

	return 1;
}








int Read_VCF_IN_Phase(In3str1v * paraFA04, Para_18 * para_18 , map <string,map <llong, vector <BaseType> > > &  SNPList ,map <string,llong> & MaxSNPSite,int & Flag_for_pro )
{
	igzstream VCFIN ((paraFA04->InStr1).c_str(),ifstream::in);

	if (VCFIN.fail())
	{
		cerr << "open VCF File IN File error: "<<(paraFA04->InStr1)<<endl;
		return  0;
	}


	vector<string> inf ;

	while(!VCFIN.eof())
	{
		string  line ;
		getline(VCFIN,line);
		if (line.length()<=0 )  { continue  ; }
		else if ( line[0] == '#' && line[1] == '#' )  { continue  ; }
		else if ( line[0] == '#' && line[1] != '#' )
		{
			inf.clear();
			split(line,inf," \t");
			if  ( inf[0]  != "#CHROM")
			{
				continue  ;
			}			
			break ;
		}
		else if ( line[0] != '#' && line[1] != '#' )
		{
			cerr<<"wrong Line : "<<line<<endl;
			cerr<<"VCF Header same thing wrong, can find sample info before site info"<<endl;
			cerr<<"VCF Header sample info Flag : [  #CHROM  ] "<<endl;
			return  0;
			break;
		}
	}

	int BadSite=0 ;
	int BadIndelSite=0;

	int Asample=inf.size();
	int SampleNum=(Asample-9);

	if (SampleNum<10)
	{
		cerr<<"sub Group Population szie is too small, recommends SampleSize >=10, SubGroup sample size: "<<SampleNum<<endl;
		return  0;
	}

	vector<string> Alt ;
	vector<string> Btmp  ;
	int  Base_len ;
	llong Site ;
	int Het_count=0;
	int Miss_count=0;

	int BaseConut=0;
	char best_base='N';
	char sed_base='N';
	int Max=0;
	int SeD=0;
	map <char,int>  :: iterator it ;
	int ERA=0;
	vector<char>  genotype ;
	vector<BaseType>  genotypeVE ;

	while(!VCFIN.eof())
	{
		string  line ;
		getline(VCFIN,line);
		if (line.length()<=0)  { continue  ; }

		istringstream isoneLine (line,istringstream::in);
		for (int iik=0 ; iik<Asample ; iik++)
		{
			isoneLine>>inf[iik];
		}
		Base_len=inf[3].length();
		Alt.clear();
		split(inf[4],Alt,",");

		for (int ii=0 ; ii<Alt.size() ;ii++)
		{
			if (Alt[ii].length()>Base_len)
			{
				Base_len=Alt[ii].length();
			}
		}

		if (Base_len>1)
		{
			BadIndelSite++;
			continue ;
		}

		map <char,int > Count ;
		Het_count=0;
		Miss_count=0;

		for (int jj=9 ; jj< Asample ;jj++ )
		{
			Btmp.clear();
			split(inf[jj], Btmp,":");
			string Genotype=Btmp[0];
			if (  Genotype[0] == '.' )
			{
				Miss_count++ ;
			}
			else
			{
				if  (Genotype[0] != Genotype[2] )
				{
					Het_count++;
				}
				Count[Genotype[0]]++;
				Count[Genotype[2]]++;
			}
		}

		if ( ( (Miss_count*1.0/SampleNum)  >(para_18->Miss)  )  )
		{
			continue ;
		}

		if ( ( (Het_count*1.0/SampleNum)  >(para_18->Het) )  )
		{
			continue ;
		}

		BaseConut=0;
		best_base='N';
		sed_base='N';
		Max=0;
		SeD=0;		

		for ( it=Count.begin(); it!=Count.end(); it++ )
		{
			if ( (it->first ) == 'N' )
			{					
				continue ;
			}
			else if ((it->second)  > Max )
			{
				SeD=Max;
				sed_base=best_base;
				Max=(it->second);
				best_base=it->first;
			}
			else if ( (it->second)  >= SeD )
			{
				SeD=(it->second);
				sed_base=it->first;
			}
			BaseConut++;
		}
		if (BaseConut==1 || BaseConut >2 )
		{
			BadSite++;
			continue ;
		}

		if ( (SeD*1.0/(SeD+Max))  < (para_18->MAF) )  
		{
			continue ;
		}

		genotype.clear();

		for (int jj=9 ; jj< Asample ;jj++ )
		{
			Btmp.clear();
			split(inf[jj], Btmp,":");
			string Genotype=Btmp[0];
			if (  Genotype[0] == '.' )
			{
				genotype.push_back('N');
				genotype.push_back('N');
			}
			else
			{
				genotype.push_back(Genotype[0]);   // phase VCF
				genotype.push_back(Genotype[2]);   // phase VCF
			}
		}



		genotypeVE.clear();

		BaseType  TypeA;

		ERA=genotype.size();

		for (int hh=0 ; hh<ERA ;hh++)
		{

			if (genotype[hh] == best_base  )
			{
				TypeA.Value=0;
			}
			else if (genotype[hh] == sed_base )
			{
				TypeA.Value=1;
			}
			else
			{
				TypeA.Value=2; 
			}
			genotypeVE.push_back(TypeA);
		}





		istringstream isone (inf[1],istringstream::in);
		isone>> Site ;


		map <string,map <llong, vector <BaseType> > >  :: iterator itSNP=SNPList.find(inf[0]);

		if (itSNP == SNPList.end())
		{
			map <llong, vector <BaseType> > DD;
			DD[Site]=genotypeVE;
			SNPList.insert(map <string,map <llong,vector <BaseType> > > ::value_type(inf[0],DD));
			Flag_for_pro++;
			MaxSNPSite.insert(map <string,llong>::value_type(inf[0],Site));
		}
		else
		{
			(itSNP->second).insert(map <llong, vector <BaseType> >  :: value_type(Site,genotypeVE)) ;
			Flag_for_pro++;
			map <string, llong> :: iterator  MaxIT =MaxSNPSite.find(inf[0]);
			if ( (MaxIT->second) < Site )
			{
				MaxIT->second=Site ;
			}

		}
	}


	VCFIN.close();

	if(BadIndelSite!=0)
	{
		cout<<"warning skip Indel site, there are total skip Indel sites number is : "<<BadIndelSite<<endl;
	}
	if (BadSite!=0)
	{
		cout<<"Warning skip non bi-allelic(Singleton/ThreeMulti allelic) site, and total skip allelic sites number is :"<<BadSite<<endl;
	}

	return 1;
}






int Identify_subgroup( string ListFile ,  map <string,int > & SubVetor  , map <string,int >  & GroupID  )
{
	igzstream SampleList (ListFile.c_str(),ifstream::in);
	if (SampleList.fail())
	{
		cerr << "open Sub Group IN File error: "<<ListFile<<endl;
		return 0;
	}

	map <string ,int >  :: iterator it;
	map <string,int > :: iterator IT_GID;

	int subpopNum=0;
	int Thisub=0;
	while(!SampleList.eof())
	{
		string  line ;
		getline(SampleList,line);
		if (line.length()<=0 || line[0] == '#' )  { continue ;}
		if (line.find(':') ==string::npos)
		{
		}
		else
		{
			vector<string> tee ;
			split(line,tee," \t:");
			line=line.substr(line.find(':')+1);
			IT_GID=GroupID.find(tee[0]);
			if ( IT_GID==GroupID.end() )
			{
				subpopNum++;
				Thisub=subpopNum ;
				GroupID.insert(map <string,int> ::value_type(tee[0],subpopNum));
			}
			else
			{
				Thisub=IT_GID->second;
			}
		}
		vector<string> inf;
		split(line,inf," \t");
		int A=inf.size();
		for (int ii=0 ; ii<A; ii++)
		{
			if (inf[ii] == ":")
			{
				cerr <<"some thing wrong for line double \':\' in lines :"<<line<<endl;
			}
			it=SubVetor.find(inf[ii]);
			if (it==SubVetor.end())
			{
				SubVetor.insert(map <string ,int> ::value_type(inf[ii],Thisub));
			}
		}
	}
	SampleList.close();

	if (subpopNum>2)
	{
		cerr <<"\tIdentifying multiple groups : ";
		for (it =GroupID.begin();it!=GroupID.end(); it++ )
		{
			cerr<<" "<<it->first;
		}
		cerr<<"\n\tThe software only accepts up to two groups. but Identify by symbol ':' Number is :" <<subpopNum<<endl;
		return   subpopNum ;
	}
	else if (subpopNum==2)
	{
		int NumA=0;
		int NumB=0;
		for (it = SubVetor.begin();it!=SubVetor.end(); it++)
		{
			if ((it->second)==1)
			{
				NumA++;
			}
			else if ((it->second)==2)
			{
				NumB++;
			}
		}

		for (it =GroupID.begin();it!=GroupID.end(); it++ )
		{
			if ((it->second)==1) 
			{
				cout << "subgroup ["<<it->first<<"] sample size is "<<NumA<<endl;
			}
			else
			{
				cout << "subgroup ["<<it->first<<"] sample size is "<<NumB<<endl;
			}
		}

		if  (NumA<10||NumB<10|| (NumA*1.0/NumB>3) || (NumB*1.0/NumA>3) )
		{
			cerr <<"In order to get the selection region by compair the two groups, it is recommended that the sample size of the two groups should not differ too much.(0.3--3) and ever subgroup sample size should >=10 ;  exit(1) \n";
			return 3;
		}
		else
		{

			return 2;
		}
	}
	else
	{
		return 1 ;
	}
}


//////////
#endif //
///////// swimming in the sky and flying in the sea ////////////


