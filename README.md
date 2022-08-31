以下の２編の論文の数値例を導出するためのMATLAB codeです。

["A design method of model error compensator for systems with polytopic-type uncertainty and disturbances(2021)"](https://www.tandfonline.com/doi/full/10.1080/18824889.2021.1918392?src=)
["Model Error Compensator Design for Continuous- and Discrete-Time Non-minimum Phase Systems with Polytopic-Type Uncertainties(2022)"](https://www.tandfonline.com/doi/pdf/10.1080/18824889.2022.2052628?needAccess=true)

pso.mで粒子群最適化を用いた設計、iterative_design.mで共通リアプノフ行列と誤差補償器の反復設計ができます。

- pso.mでは、kMaxが最大更新回数、numParticlesが粒子数、omega,c1,c2が更新式の重みになります。

- iterative_design.mを実行する場合、初期値としてAd,Bd,Cd,Ddの値が必要になります。initialValueFilePathでresultsの中のファイルを指定すれば、PSOの設計結果を用いることもできますが、Ad,Bd,Cd,Ddを手動で与えることもできます。

- plants/stable_mp.matに最小位相系の制御対象の、plants/stable_nmp.matに非最小位相系の制御対象の情報が保存されています。制御対象を変えたい場合はこれらのファイルの中身の行列の値を変更します。(A,B,C)が制御対象、(Am,Bm,Cm)がモデル、(Af,Bf,Cf)がPFCの動特性を表す行列です。  
状態数をn、入力数をi、出力数をo、端点数をvとすると、Aはn\*n\*v、Bはn\*i\*v、Cはo\*n\*v、Bwはn\*i、Dwはo\*oの行列になります。

- 粒子群最適化で設計する際に誤差補償器の構造を変更したい場合は、rAd,rBd,rCd,rDdを変更します。  
反復設計で設計する際に誤差補償器の構造を変更したい場合は、search_d.mもしくはsearch_d_pfc.mのcanonical_formの値を変更します。また、可制御正準形で次数やその構造を変更したい場合は、関数search_d.mもしくはsearch_d_pfc.mの中のTal,Tar,Tas,Tbr,Tbsの構造とsetlmisの下の\[Ad,~,sAd\] = lmivar(2,\[2,2]);と、\[Bd,\~,sBd\] = lmivar(2,\[2,2]);を変更します。

- 設計した誤差補償器を用いた応答を確認する場合は、test_response.mを実行します(simulinkが必要です)。設計結果がresults/にあるので、matFileでそのファイルのパスを指定します。また、PFCを使って設計した場合はusePFC=trueにします。

以下各ファイルの概要です。

- plants/  
    MECを設計する制御対象をまとめたディレクトリ。  
    - stable_mp.mat  
        最小位相系の制御対象。["A design method of model error compensator for systems with polytopic-type uncertainty and disturbances(2021)"](https://www.tandfonline.com/doi/full/10.1080/18824889.2021.1918392?src=)の数値例で用いた制御対象。
    - stable_nmp.mat  
        非最小位相系の制御対象。["Model Error Compensator Design for Continuous- and Discrete-Time Non-minimum Phase Systems with Polytopic-Type Uncertainties(2022)"](https://www.tandfonline.com/doi/pdf/10.1080/18824889.2022.2052628?needAccess=true)の数値例で用いた制御対象。
- results/  
    iteretave_design.m、pso.mの設計結果を保存するディレクトリ。  
- simulink_models/  
    応答の確認に用いるsimulinkのモデルが入ったディレクトリ。
- iterative_design.m  
    誤差補償器と共通リアプノフ行列の反復設計を行うスクリプト。
- pso.m  
    粒子群最適化で誤差補償器を設計するスクリプト。
- search_d.m  
    反復設計で用いる、共通リアプノフ行列を与えてガンマを最小化する誤差補償器のパラメータを得る関数(PFCなし)。
- search_x.m  
    反復設計で用いる、誤差補償器を与えて最小のガンマとその時の共通リアプノフ行列の値を得る関数(PFCなし)。
- search_d_pfc.m  
    反復設計で用いる、共通リアプノフ行列を与えてガンマを最小化する誤差補償器のパラメータを得る関数(PFCあり)。
- search_x_pfc.m  
    反復設計で用いる、誤差補償器を与えて最小のガンマとその時の共通リアプノフ行列の値を得る関数(PFCあり)。
- test_response.m  
    粒子群最適化もしくは反復設計で設計した誤差補償器を制御対象に適用したときの応答を確認するスクリプト。
