# PlintoPali - Analisi Statica e Sismica di Plinti su Pali

Applicazione web interattiva basata su Streamlit per il calcolo e la verifica di fondazioni profonde (plinti su pali). Il tool esegue l'analisi della capacità portante del palo singolo e calcola la ripartizione dei carichi sui pali del gruppo in condizioni statiche e sismiche, tenendo conto della stratigrafia del terreno e della presenza di falda acquifera.

## 📖 Fondamenti Teorici

Il motore di calcolo dell'applicazione si basa su principi consolidati di geotecnica e meccanica delle strutture.

### 1. Tensioni Efficaci e Falda
Le tensioni verticali efficaci $\sigma'_{v}$ vengono calcolate integrando i pesi di volume del terreno. In presenza di falda acquifera, il peso dell'unità di volume viene decurtato della spinta idrostatica, utilizzando il peso di volume alleggerito $\gamma' = \gamma_{sat} - \gamma_{w}$.

### 2. Capacità Portante del Palo Singolo
La capacità portante ultima ($Q_{ult}$) è calcolata come somma della portata di punta ($Q_b$) e della resistenza laterale ($Q_s$):
$$Q_{ult} = Q_b + Q_s$$

* **Resistenza di Punta ($Q_b$):** Dipende dalla natura del terreno alla base del palo.
    * Terreni a grana grossa (drenati): $Q_b = N_q \sigma'_{v}(L) A_{base}$.
    * Terreni a grana fina (non drenati): $Q_b = N_c c_u A_{base}$.
* **Resistenza Laterale ($Q_s$):** Calcolata integrando l'attrito lungo il fusto per ogni strato.
    * Metodo $\beta$ (terreni attritivi): $q_s = \beta \sigma'_{v}(z)$.
    * Metodo $\alpha$ (terreni coesivi): $q_s = \alpha c_u$.

### 3. Ripartizione dei Carichi (Plinto Rigido)
Il calcolo delle reazioni sui singoli pali assume che la piastra di fondazione (plinto) sia infinitamente rigida. La distribuzione delle sollecitazioni ($N$, $M_x$, $M_y$) avviene risolvendo il seguente sistema lineare sovradeterminato ai minimi quadrati:
$$\mathbf{R} = \mathbf{A} (\mathbf{A}^T \mathbf{A})^{-1} \mathbf{F}$$
Dove $\mathbf{A}$ è la matrice delle coordinate geometriche dei pali, $\mathbf{F}$ è il vettore delle azioni (modificate per includere le componenti pseudo-statiche sismiche $k_h$ e $k_v$) e $\mathbf{R}$ è il vettore delle reazioni verticali sui pali.

## ⚙️ Dati di Input e Output Attesi

### Input
L'applicazione accetta i seguenti parametri tramite barra laterale:
* **Geometria Plinto:** Dimensioni $B$ ed $L$, numero di pali per direzione ($n_x$, $n_y$), e interassi.
* **Geometria Palo:** Diametro e lunghezza.
* **Azioni di Progetto:** Sforzo normale $N$ e momenti flettenti $M_x$, $M_y$.
* **Parametri Geotecnici:** Coefficienti di portanza ($N_q$, $N_c$), coefficienti di attrito ($\alpha$, $\beta$) e fattore di sicurezza globale.
* **Azione Sismica e Falda:** Coefficienti sismici $k_h$ e $k_v$, e profondità della falda dal piano campagna.
* **Stratigrafia:** Formato CSV multiriga contenente spessore, $\gamma_{dry}$, $\gamma_{sat}$, angolo di attrito $\phi$, coesione non drenata $c_u$ e modulo $k$.

### Output
* **Metriche di Sintesi:** Numero totale dei pali, portata ammissibile, e fattori di sicurezza minimi.
* **Avvisi e Note Tecniche:** Controlli automatici su trazione dei pali, superamento della portata e interferenze da "effetto gruppo".
* **Tabelle:** Sintesi dei parametri e tabella comparativa delle reazioni statiche vs sismiche esportabile in CSV.
* **Visualizzazioni Plotly:** Pianta del plinto con heatmap delle reazioni, diagrammi a bolle per i carichi e grafici a barre comparativi.

## ✨ Caratteristiche Principali

* **Interfaccia Intuitiva:** Sviluppata interamente in Streamlit per un utilizzo immediato tramite browser.
* **Import/Export:** Salvataggio e ricaricamento degli input in formato JSON.
* **Analisi Comparativa:** Valutazione simultanea delle combinazioni statiche e pseudo-statiche (sismiche).
* **Grafica Avanzata:** Utilizzo di `plotly.graph_objects` per rendering di geometrie e distribuzioni di carico.

## 🛠 Requisiti di Sistema

* Python 3.8 o superiore
* Librerie: `streamlit`, `numpy`, `pandas`, `plotly`

## 🚀 Installazione

1.  Clona il repository:
    ```bash
    git clone https://github.com/tuouser/plintopali.git
    cd plintopali
    ```
2.  Installa le dipendenze:
    ```bash
    pip install -r requirements.txt
    ```

## 💻 Utilizzo

Per lanciare l'applicazione localmente, esegui il seguente comando dalla radice del progetto:
```bash
streamlit run app.py
```
L'interfaccia si aprirà automaticamente nel tuo browser predefinito, solitamente all'indirizzo `http://localhost:8501`.

## 🤝 Contributi

I contributi sono benvenuti! Sentiti libero di aprire una Issue per segnalare bug o proporre nuove funzionalità, oppure invia una Pull Request con le tue migliorie.
