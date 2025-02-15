#include<bits/stdc++.h>
using namespace std;
int pop=100;
const int T=100;
int N;
int NE;
const double p=1.0;
const double lenP=5.0;
int Ne=2;
vector<vector<bool>> A;
vector<vector<int>> x(pop+1);
vector<vector<int>> e;
vector<int> xBest;
vector<int> d;
vector<vector<int>> dk(pop+1);
vector<vector<int>> lk(pop+1);
vector<int> pos(Ne);
random_device rd;   
mt19937 gen(rd());
double rateInc=0.1;
double rateDec=0.3;
double rateMimicElist=0.8;
int maLong=0;
int macom=0;
double calI(vector<int> l1,vector<int> l2){
    int s1 = *max_element(l1.begin(), l1.end());
    int s2 = *max_element(l2.begin(), l2.end());
    vector<vector<int>> c1 (s1+1);
    vector<vector<int>> c2(s2+1);

    for (int i=1;i<=N;i++)
        c1[l1[i]].push_back(i),c2[l2[i]].push_back(i);

    double I=0;
    for (int i=1;i<=s1;i++){
        if (c1[i].empty()) continue;
        for (int j=1;j<=s2;j++){
            if (c2[j].empty()) continue;
            
            double Cij=0;
            for (int v:c1[i]) if (l2[v]==j) Cij++;
            if (Cij>0)
                I+=Cij*log( Cij*double(N)/double(c1[i].size()*c2[j].size()) );
        }
    }

    return I;
}
double calH(vector<int> l){
    int s=*max_element(l.begin(), l.end());
    vector<vector<int>> c(s+1);
    for (int i=1;i<=N;i++)
        c[l[i]].push_back(i);
    double H=0;
    for (int i=1;i<=s;i++){
        if (c[i].size())
            H+=double(c[i].size()*log(double(c[i].size())/double(N)));
    }
    return H;
}
double NMI(vector<int> l1,vector<int> l2){
    return -2.0*calI(l1,l2)/(calH(l1)+calH(l2));
}

double modularity(vector<int> dk,vector<int> lk){
    double Q=0;

    for (int i=1;i<=N;i++){
        Q+=double(lk[i])/double(NE)-pow(double(dk[i])/double(2*NE),2.0);
    }

    return Q;
}

void transfer(vector<int> &dk,vector<int> &lk,vector<int> l,int i,int l1,int l2){
    dk[l1]-=d[i];
    dk[l2]+=d[i];
    for (int v:e[i]){
        if (l1==l[v])
            --lk[l1];
        if (l2==l[v])
            ++lk[l2];
    }
}
void LAR_rand(vector<vector<int>> &a){
    for (int u=1;u<=N;u++){
        if (!e[u].size()) continue;
        uniform_int_distribution<int> disv(0,e[u].size()-1);
        int v=e[u][disv(gen)];
        a[u].push_back(v);
        a[v].push_back(u);
    }
}

vector<int> decoding(vector<vector<int>> a){
    bool dd[N+1]={};
    vector<int> l(N+1);
    int cnt=0;

    for (int i=1;i<=N;i++)
        if (!dd[i]){
            ++cnt;
            queue<int> q;
            q.push(i);
            while (!q.empty()){
                int u=q.front();
                q.pop();
                l[u]=cnt;
                for (int v:a[u])
                    if (!dd[v]){
                        dd[v]=true;
                        q.push(v);
                    }
            }
        }
    return l;
}

void initialization(){
    int s;
    for (int p=1;p<=pop;p++){
        vector<vector<int>> a(N+1);
        LAR_rand(a);
        x[p]=decoding(a);         
        
        s=*max_element(x[p].begin(),x[p].end());
        dk[p].resize(N+1,0);
        lk[p].resize(N+1,0);                            
        
        for (int u=1;u<=N;u++){
            dk[p][x[p][u]]+=d[u];
             
            for (int v:e[u])
                if (x[p][u]==x[p][v]&&u<v){
                    ++lk[p][x[p][u]];
                }
        }
    }    
}

void movingToPrey(vector<int> &l,vector<int> &dk,vector<int> &lk,double k,vector<int> xHybrid){
    vector<int> pos;
    for (int i=1;i<=N;i++)
            pos.push_back(i);

    shuffle(pos.begin(),pos.end(),gen);

    for (int i=0;i<=int(k)-1;i++){
        transfer(dk,lk,l,pos[i],l[pos[i]],xBest[pos[i]]);
        l[pos[i]]=xBest[pos[i]];
    }
}

void randomWalk(vector<int> &l,vector<int> &dk,vector<int> &lk,double k){
    vector<int> ranPop;
    for (int i=1;i<=pop;i++)
        ranPop.push_back(i);
    shuffle(ranPop.begin(),ranPop.end(),gen);   

    vector<vector<int>> P;
    for (int i=0;i<=lenP-1;i++)
        P.push_back(x[ranPop[i]]);
    

    vector<int> randNode;
    for (int i=1;i<=N;i++)
        randNode.push_back(i);
    shuffle(randNode.begin(),randNode.end(),gen);


    for (int i=lenP-1;i>=0;i--){
        int kCommma=int(k)/(i+1)+int(int(k)%(i+1)>0);

        while (kCommma--){
            --k;
            transfer(dk,lk,l,randNode[k],l[randNode[k]],P[i][randNode[k]]);
            l[randNode[k]]=P[i][randNode[k]];
        }
    }
}

void encirlingThePrey(vector<int>&l,vector<int> &dk,vector<int> &lk,double r,vector<int> xHybrid){
    int k=int(r*double(N));

    vector<int> pos;
    for (int i=1;i<=N;i++)
            pos.push_back(i);

    shuffle(pos.begin(),pos.end(),gen);

    for (int i=0;i<=int(k)-1;i++){
        transfer(dk,lk,l,pos[i],l[pos[i]],xBest[pos[i]]);
        l[pos[i]]=xBest[pos[i]];
    }
}

void mutation(vector<int> &l,vector<int> &dk,vector<int> &lk,double u){
    uniform_real_distribution<double> dis(0,1);
    
    int S = *max_element(l.begin(), l.end());
    vector<int> ltmp;
    vector<int> dktmp;
    vector<int> lktmp;
    for (int i=1;i<=N;i++){
        double x=dis(gen);
        if (x<u){
            ltmp=l; 
            dktmp=dk;
            lktmp=lk;

            ++S;
            double y=dis(gen);
            if (y<0.5){
                transfer(dk,lk,l,i,l[i],S);                
                l[i]=S;
                
            }
            else{
                transfer(dk,lk,l,i,l[i],S);
                l[i]=S;
                

                for (int neigbor:e[i]){
                    transfer(dk,lk,l,neigbor,l[neigbor],S);
                    l[neigbor]=S;
                }
            }
            
            if (modularity(dk,lk)<modularity(dktmp,lktmp)){
                --S;
                l=ltmp;
                dk=dktmp;
                lk=lktmp;
            }
            macom=max(macom,S);
        }
    }
}


void boudaryNodeAdjustment(vector<int> &l,vector<int> &dk,vector<int> &lk,double rate){
    vector<int> tmpl;
    vector<int> dktmp;
    vector<int> lktmp;
    int s=*max_element(l.begin(),l.end());
    bool dd[s+1]={};
    
    uniform_real_distribution dis(0.0,1.0);
    for (int i=1;i<=N;i++){
            if (dis(gen)>rate) continue;

            dd[l[i]]=true;
            for (int neighbor:e[i])
                if (l[i]!=l[neighbor] && !dd[l[neighbor]]){
                    dd[l[neighbor]]=true;
                    
                    tmpl=l;
                    dktmp=dk;
                    lktmp=lk;

                    transfer(dk,lk,l,i,l[i],l[neighbor]);
                    l[i]=l[neighbor];

                    if (modularity(dk,lk)<modularity(dktmp,lktmp)){
                        dk=dktmp;
                        lk=lktmp;
                        l=tmpl;
                    }
                }

            for (int i=1;i<=s;i++)
                dd[i]=false;
        }
    
}
void EPD(){
    if (x.size()<10) return;

    vector<pair<double, int>> modularityValues;
    for (int i = 1; i <= pop; i++) {
        double modValue = modularity(dk[i],lk[i]);  
        modularityValues.push_back({modValue, i});
    }

    sort(modularityValues.begin(), modularityValues.end());
    reverse(modularityValues.begin(),modularityValues.end());

    vector<vector<int>> sortedX(pop + 1);
    vector<vector<int>> sorteddk(pop + 1);
    vector<vector<int>> sortedlk(pop + 1); 
    for (int i = 0; i < pop; i++) {
        sortedX[i + 1] = x[modularityValues[i].second];
        sorteddk[i + 1] = dk[modularityValues[i].second];
        sortedlk[i + 1] = lk[modularityValues[i].second];
    }

    x=sortedX;
    dk=sorteddk;
    lk=sortedlk;
    
    // find EL
    vector<vector<int>> EB;
    for (int i=1;i<=Ne;i++) EB.push_back(x[i]);

    int cnt=0;

    for (vector<int> xb:EB){
        double suit,miSuit=1;
        for (int i=3;i<=pop;i++){
            suit=NMI(x[i],xb);
            if (suit<miSuit){
                miSuit=suit;
                pos[cnt]=i;
            }

            if (suit==1){
                for (int j=0;j<cnt;j++)
                    pos[j]-=(i<pos[j]);

                x.erase(x.begin() + i);
                dk.erase(dk.begin() + i);
                lk.erase(lk.begin() + i);   
                --pop;
                --i;
            }
        }
        ++cnt;
    }

    //real EPD
    double N_nor=pop-(pop/2+1)+1;
    uniform_real_distribution<double> dis(0,1);
    for (int i=pop/2+1;i<=pop;i++)
        if (i!=pos[0]&&i!=pos[1]){
            double C=1.0-exp(-double(i)/N_nor);
            double rand=dis(gen);
            if (rand<=C){
                for (int j=0;j<Ne;j++)
                    pos[j]-=(i<pos[j]);

                x.erase(x.begin() + i);
                dk.erase(dk.begin() + i);
                lk.erase(lk.begin() + i);   
                --pop;
                --i;
            }
        }
    if (Ne)
        maLong=NMI(x[pos[Ne-1]],xBest);
}

void updateLocation(vector<int> &l,int t,vector<int> &dk,vector<int> &lk){
    uniform_real_distribution<double> dis(0,1);
    double alpha=dis(gen),beta=dis(gen);
    if (alpha<0.5){
        double k=p*double(t)*double(N)/double(T);

        if (beta<0.5){
            if (dis(gen)<rateMimicElist)
                movingToPrey(l,dk,lk,k,xBest);
            else movingToPrey(l,dk,lk,k,x[pos[0]]);

        }
        else
            randomWalk(l,dk,lk,k);
    }
    else {
        uniform_real_distribution<double> disl(-1,1);
        double ll=disl(gen);
        double r=abs(cos(2*M_PI*ll));

        if (dis(gen)<rateMimicElist)
            encirlingThePrey(l,dk,lk,r,xBest);
        else  encirlingThePrey(l,dk,lk,r,x[pos[0]]);
    }
}

void EP_WOCD(){
    initialization();
    double ans=0;
    for (int i=1;i<=pop;i++)
        if (modularity(dk[i],lk[i])>ans){
            ans=modularity(dk[i],lk[i]);
            xBest=x[i];
        }
    
    uniform_real_distribution dis(0.0,1.0);
    for (int t=1;t<=T;t++){
        for (int p=1;p<=pop;p++){
            double rateLS=1,rateMu=0.3;
            bool check=(p>Ne);
            for (int i=0;i<Ne;i++)
                if (p==pos[i]) {check=0;rateLS=rateMu=0.1;break;}

            if (check) updateLocation(x[p],t,dk[p],lk[p]);
            mutation(x[p],dk[p],lk[p],rateMu);
            boudaryNodeAdjustment(x[p],dk[p],lk[p],rateLS);
        }

        bool isStable=1;
        for (int i=1;i<=pop;i++){
            if (modularity(dk[i],lk[i])>ans){
                ans=modularity(dk[i],lk[i]);
                xBest=x[i];

                isStable=0;
            }
        }

        // if (!isStable){
        //     if (dis(gen)<rateInc){
        //         --Ne;
        //     }
        // }
        // else {
        //     if (dis(gen)<rateDec){
        //         ++Ne;
        //     }
        // }

        // Ne=min(Ne,pop/2);

        cout<<ans<<" "<<Ne<<" "<<maLong<<" "<<macom<<"\n";
        EPD();        
    }    

    cout<<ans<<"\n";
    for (int i=1;i<=N;i++)
        cout<<xBest[i]<<" ";
    cout<<"\n";
    cout<<pos[Ne-1]<<" "<<pos[Ne-2];
}
int main(){
    clock_t tStart = clock();

    freopen("input.txt","r",stdin);

    cin>>N;
    cin>>NE;

    d.resize(N+1);
    A.resize(N + 1, vector<bool>(N + 1, 0));
    e.resize(N+1);
    for (int i=1;i<=NE;i++){
        int u,v;
        cin>>u>>v;
        u++,v++;
        e[u].push_back(v);
        e[v].push_back(u);
        d[u]++,d[v]++;
        A[u][v]=A[v][u]=true;
    }

    EP_WOCD();

    printf("\nTime taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
}