#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include "MurmurHash3.h"
#include <string.h>
#include <math.h>
#include <vector>
#include <cstdint>
#include <limits>
#include <map>
#include<bits/stdc++.h>
#include <utility>
#include <type_traits>
#include <typeinfo>
#include <cxxabi.h>
#include <chrono>

using namespace std;
using namespace std::chrono;

int main(){
    ifstream inFile;
    inFile.open("ERR626208.fastq");
    string line;
    //const int d=4;  // filas de funciones hash.
    //const int w=65536;  // # de buckets para cada fila. ///{2^15= 32768 y 2^16=65536}
    const int q=20;     ///{q-mer size: 20, 30, 40}
    const int limit=166-q;
    //const int UMBRAL=16000; // umbral de frecuencia para HH´s. //{CountminCU entre 5000 y 8000, Countmin 16000}
    //const int DEFAULT=0;
    std::vector<string> elem;
    std::vector<int> cont;
    std::vector<string>::iterator it;
    int k=250;  //  1/phi
    //unordered_map<string, int> hhitters;
    //unordered_map<string, int>:: iterator itr;
    int cont_line=0;
    bool flag=false;
    //int cont_N;
    // 1era pasada para calcular HH's
    auto start = high_resolution_clock::now();
    while (inFile){ // llenado de matriz de contadores.
        getline(inFile, line);
        if(line[0]=='@'){
            continue;
        }
        else if(line[0]=='+'){
            flag=true;
            continue;
        }
        else if(flag==true){
            flag=false;
            continue;
        }
        else{
            cont_line++;
            if(line.size()==0){
                break;
            }
            for(int i = 0; i <= limit; i++){
                string temp = line.substr(i, q);   //  substr(i, q-mer size) ///{q-mer size: 20, 30, 40}
                //char qmer[temp.size()+1];
                //strcpy(qmer, temp.c_str());
                    //std::cout << "fafa1" << '\n';
                    it = std::find(elem.begin(), elem.end(), temp);
                    if(it != elem.end()){
                        //std::cout << "fafa2" << '\n';
            			int index = it - elem.begin();
            			cont[index]++;
            		}
                    else if(elem.size() < k){
                        //std::cout << "fafa3" << '\n';
                        elem.push_back(temp);
                        cont.push_back(1);
                    }
                    else{
                        //std::cout << "fafa4" << '\n';
                        for(int l = 0; l < k; l++){
                            cont[l]--;
                        }
                    }
                    for (int l = 0; l < elem.size(); l++) {
                        if(cont[l] <= 0){
                            cont.erase(cont.begin() + l);
                            elem.erase(elem.begin() + l);
                        }
                    }

            }
       }//std::cout << cont_line << '\n';
    }
    //cont_N = (166-q+1)*cont;
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    std::cout << '\n' << "Ejecucion: "<< (duration.count())/60 << " minutos, "<< (duration.count())%60 << " segundos.\n";
    cout << "\nLineas procesadas: "<< cont_line << endl;
    //std::cout << "N= " << cont_N <<'\n';
    cout << "\n" << "Resultados con MisraGries, " << "Parametros: k = "<< k << ", q = "<< q <<"\n";
    std::cout << "\nTamaño del vector de elementos: " << elem.size() << '\n' << endl;
    for (int i = 0; i < elem.size(); i++) {
        std::cout << elem[i] << " " << cont[i] << '\n';
    }
    //cout << "\n" << "Resultados con misraGries: " << hhitters.size() << " H.Hitters con umbral de " << UMBRAL << "\n";

    /*for (itr = hhitters.begin(); itr != hhitters.end(); itr++){
        cout << itr->first << "  " << itr->second << endl;
    }
    /*inFile.clear();
    inFile.seekg(0,std::ios::beg);
    // 2da pasada
    int freq_est;
    while (inFile){
        getline(inFile, line);
        if(line[0]=='@'){
            continue;
        }
        else if(line[0]=='+'){
            flag=true;
            continue;
        }
        else if (flag==true){
            flag=false;
            continue;
        }
        else{
            if(line.size()==0){
                break;
            }
            for(int i = 0; i <= 146; i++){
                string temp = line.substr(i, 20);   //  substr(i, q-mer size)
                char qmer[temp.size()+1];
                strcpy(qmer, temp.c_str());
                for(int seed = 0; seed < d; seed++) {
                    MurmurHash3_x86_32(qmer, (uint16_t)strlen(qmer), seed, hash_otpt);//valor hashed
                    if(seed == 0){   // se busca el minimo en las d filas para el valor hashed
                        freq_est = C[seed][hash_otpt[0]];
                    }
                    else{
                        if(C[seed][hash_otpt[0]] < freq_est){
                            freq_est = C[seed][hash_otpt[0]];
                        }
                    }
                }
                if(freq_est >= UMBRAL){ //  agregar HH´s por sobre UMBRAL = 15000
                    hhitters[qmer] = freq_est;
                }
            }
       }
    }
    /*for (int d_i = 0 ; d_i < C.size () ; d_i++) {
        for (int w_i = 0 ; w_i < C[d_i].size () ; w_i++) {
            if(C[d_i][w_i]>=10000){
                cout << C[d_i][w_i] << ' ' ;
            }
        }
    }*/
    inFile.close();
    return 0;
}
