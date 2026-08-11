# L'instabilità essenziale del vuoto di Yang–Mills: riproduzione machine-checked del programma di Preparata

Versione 1.1 — 11 agosto 2026 (v1.1: testo primario di Preparata acquisito e verificato a livello di equazioni — Sezione 8)

Giuliano (Vantasner AG), per Luca Gamberale

Rapporto tecnico con verifica automatica (substrate-framework, tracker #44)

## Sommario

Questo rapporto riproduce, con oracoli macchina (SymPy/SciPy/mpmath, 63 test automatici), il programma di Giuliano Preparata sull'instabilità del vuoto perturbativo di Yang–Mills: il potenziale efficace a un loop nel fondo cromomagnetico costante (vuoto di Savvidy), il potenziale a due loop in funzione di $b/\Lambda^2$ (con $b = gH$), l'analisi di stabilità rispetto alle fluttuazioni, la decomposizione per colori e per modi trasversi, una configurazione classica migliore di quella di Savvidy (il vuoto "spaghetti" di Ambjørn–Olesen) con la riquantizzazione attorno ad essa, e la verifica equazione per equazione del testo primario del 1986 (Sezione 8). Ogni affermazione quantitativa è prodotta da codice eseguibile nel repository `vantasnerdan/substrate-framework` (issue #44–#50 e #56, PR #51–#55 e #57) o citata esplicitamente dalla letteratura. Risultato principale: il risultato a un loop di Savvidy si riproduce esattamente; il risultato a due loop pubblicato (Bordag–Skalozub 2022) contiene tre errori tipografici/sostanziali che documentiamo con valori esatti; la configurazione di Savvidy è instabile alle fluttuazioni a ogni ordine perturbativo accessibile, e il condensato di tubi di flusso (reticolo triangolare, parametro di Abrikosov $\beta = 1.159595$) è energeticamente preferito.

## 1. Fonti e accessibilità

Dichiariamo esattamente cosa è stato verificato e da dove.

- G. Preparata, *Essential Quantum Instability of the Perturbative Yang–Mills Vacuum*, Nuovo Cim. A **96** (1986) 366, DOI 10.1007/BF02833896 — **in nostro possesso** (PDF di 28 pp. con testo estraibile, md5 34f8efc4c056837869cfbe63906ea14d), acquisito l'11 agosto 2026 e fissato come provenance della campagna in `proposals/P229-preparata-qcd-vacuum-audit/sources/`. La scansione del 1986 degrada le formule in esposizione: nessuna equazione è trascritta dall'OCR come verità; ogni affermazione verificabile è codificata dal suo enunciato testuale più chiaro e poi controllata con oracoli (Sezione 8).
- G.K. Savvidy, Phys. Lett. B **71** (1977) 133 — il risultato a un loop, ripreso dalle sue riesposizioni aperte (EPJC **80** (2020) 165, arXiv:1910.00654; la rassegna EPJC 2026). Verificato interamente.
- N.K. Nielsen, P. Olesen, Nucl. Phys. B **144** (1978) 376 — paywalled; il modo instabile e la parte immaginaria sono stati verificati indipendentemente dal nostro calcolo (zeta dell'operatore e integrale di modo).
- M. Bordag, V. Skalozub, EPJC **82** (2022) 390, arXiv:2112.01043 — **aperto**; è la fonte del due loop a $T=0$ (eq. 57–58). Verificato; tre discrepanze documentate (Sezione 5).
- M. Bordag, *Symmetry* **15** (2023) 1137 — rassegna aperta; fonte della costruzione del reticolo di tubi di flusso (eq. 146–152).
- P. Cea, arXiv:2311.14791 — aperto; settori tachiometrici $SU(3)$.

Metodo: nessun passaggio "a mano". Ogni coefficiente è prodotto da codice; ogni check ha una mutazione che lo rompe. I moduli: `chromomagnetic_background.py` (21 test), `chromomagnetic_two_loop.py` (13 test), `chromomagnetic_sectors.py` (10 test), `spaghetti_vacuum.py` (8 test), `preparata_1986.py` (11 test) — 63 in totale.

## 2. Il potenziale a un loop (Savvidy), verificato

Convenzioni congelate: segnatura mostly-plus, gauge di fondo di Feynman ($\xi = 1$), $SU(2)$ con fondo in una direzione di Cartan ($a = 3$), variabile covariante $b = gH$.

**Spettro delle fluttuazioni** (derivato dall'operatore discretizzato, non asserito): vettore carico
$$E^2 = p_z^2 + b(2n+1) - 2b\,s_3, \qquad s_3 \in \{+1, 0, 0, -1\},$$
fantasma carico $E^2 = p_z^2 + b(2n+1)$. Il livello $n=0$, $s_3 = +1$ dà $E^2 = p_z^2 - b$: il **tachione di Nielsen–Olesen**.

**Cancellazione fantasma**: la traccia di spin $2\cosh(2bs) + 2$ meno il fantasma ($-2$) lascia $2\cosh(2bs)$; il coefficiente di Seeley–DeWitt è esatto:
$$\frac{b}{\sinh(bs)}\,2\cosh(2bs) = \frac{2}{s} + \frac{11}{3}b^2 s + \frac{127}{180}b^4 s^3 + \cdots$$

**Potenziale a un loop** (MS-bar):
$$V_1(b) = \frac{11}{48\pi^2}\,b^2\left(\ln\frac{b}{\mu^2} - \frac{1}{2}\right).$$
Il coefficiente $11/(48\pi^2)$ è fissato a $10^{-10}$ dalla rotta zeta dell'operatore e a $10^{-2}$ dalla rotta proper-time con cut-off, indipendenti. Il minimo:
$$\ln\frac{b_{\min}}{\mu^2} = -\frac{24\pi^2}{11 g^2}, \qquad V_{\min} = -\frac{11\mu^4}{96\pi^2}\,e^{-48\pi^2/11g^2} < 0,$$
identico alla formula pubblicata da Savvidy. Con il running a un loop ($b_0 = 22/3$): $b_{\min} = \Lambda^2$ **esattamente** — contenuto indipendente dallo schema.

## 3. Stabilità rispetto alle fluttuazioni

Il verdetto è **negativo** e lo enunciamo con precisione:

- Per ogni $b > 0$ esiste il modo tachiometrico $E^2 = p_z^2 - b < 0$ per $p_z^2 < b$.
- Parte immaginaria a un loop, derivata dall'integrale di modo e indipendentemente dalla zeta dell'operatore:
$$|\,\mathrm{Im}\,V_1| = \frac{b^2}{8\pi}.$$
- A due loop la parte immaginaria persiste: dal quadrato complesso $(3g^2/2)B_2^2$ con $B_2 = -b(\ln 2 - i\pi)/(16\pi^2)$ si ottiene $\mathrm{Im}\,V_2 = -3g^2 b^2 \ln 2/(256\pi^3) \neq 0$.
- La risommazione ad anelli ("daisy") che rimuove la parte immaginaria a un loop è insufficiente a due loop (Bordag–Skalozub, conclusioni). Non esiste quindi un verdetto di stabilità perturbativo a questo ordine: il vuoto di Savvidy **non è stabile** contro le fluttuazioni, in accordo con la tesi di "instabilità essenziale" di Preparata.

## 4. Decomposizione per colori e modi trasversi

Per settore (coefficienti esatti del nucleo di calore, termine $b^2 s$):

| settore | coefficiente | nota |
| --- | --- | --- |
| trasverso $s_3 = +1$ | $11/6$ | contiene il tachione |
| trasverso $s_3 = -1$ | $11/6$ | |
| coppia longitudinale $s_3 = 0$ | $-1/3$ | |
| fantasma | $+1/3$ | cancella la coppia longitudinale |

L'intero logaritmo a un loop è trasverso-paramagnetico: i due settori di spin trasverso portano ciascuno $11/6$, sommando all'$11/3$ totale; settore longitudinale e fantasma si cancellano esattamente.

**$SU(3)$** (fondo lungo $\lambda_3$): le radici cariche hanno cariche $(1, \tfrac12, \tfrac12)$ in unità di $gH$ — i tre settori tachiometrici di Cea. Ogni radice contribuisce come un vettore carico complesso di $SU(2)$ con la sua carica; il coefficiente del logaritmo totale è $\tfrac32 C$ con $C = 11/(48\pi^2)$, più gli spostamenti $\ln q$ per radice.

## 5. Il due loop: verifica e tre errata documentate

Dal calcolo delle definizioni del paper (loro eq. (3)) a $T=0$, $\xi = 1$: il contributo a due loop è $W_2 = \tfrac32 g^2 B_2(0,b)^2$ con il "tadpole" magnetico
$$B_2(0,b) = \frac{b}{4\pi}\sum_{n,\sigma=\pm2}\int\frac{d^2k}{(2\pi)^2}\frac{1}{k^2 + b(2n+1+\sigma) - i0} = -\frac{b\,(\ln 2 - i\pi)}{16\pi^2},$$
esatto via funzione zeta: il polo e i termini $\ln b$, $\ln\mu$ decadono identicamente perché $G(1) = 1 + \zeta(0,\tfrac12) + \zeta(0,\tfrac32) = 0$ (la cancellazione coinvolge il tachione), e la somma residua $\sum \ln(2n+1+\sigma) = \ln 2 - i\pi$ segue dall'identità di Lerch $\zeta'(0,a) = \ln\Gamma(a) - \ln\sqrt{2\pi}$ (verificata numericamente a $10^{-30}$).

**Errata candidati** (ciascuno fissato da un test):

- **F1** — la loro eq. (58) stampa la correzione all'esponente $3\ln^2 2/(98\pi^2)$ nella forma esponenziale e $3\ln^2 2/(88\pi^2)$ nell'espansione. La minimizzazione della loro stessa (57) dà $88$: il $98$ è un refuso.
- **F2** — la loro (57) stampa la parte immaginaria a un loop $-i b^2/(8\pi^2)$; due rotte indipendenti (e Nielsen–Olesen) danno $b^2/(8\pi)$.
- **F3** — la loro (57) stampa il due loop $g^2\ln^2 2\, b^2/(128\pi^4)$. Dalle loro stesse definizioni:
$$\mathrm{Re}\,W_2 = \frac{3g^2 b^2 (\ln^2 2 - \pi^2)}{512\pi^4},$$
cioè un fattore $4/3$ sul coefficiente di $\ln^2 2$ e una struttura $-\pi^2$ assente nella forma stampata (differenza esatta: $-g^2b^2(\ln^2 2 + 3\pi^2)/(512\pi^4)$). Riprodotto indipendentemente dall'estrattore bibliografico della campagna.

**Minimo a due loop** (con il nostro $W_2$): $b_{\min} = \mu^2 \exp\left(-\frac{24\pi^2}{11g^2} - \frac{9(\ln^2 2 - \pi^2)g^2}{352\pi^2}\right)$.

## 6. Il potenziale in funzione di $b/\Lambda^2$

Con running a un loop e miglioramento RG ($\mu^2 = b$), $x \equiv b/\Lambda^2$:
$$\frac{W(x)}{\Lambda^4} = \frac{11}{48\pi^2}\,x^2\left(\ln x - \frac12\right) + \kappa\,\frac{x^2}{\ln x},\qquad \kappa = \frac{9(\ln^2 2 - \pi^2)}{704\pi^2}\ \text{(derivato)},\quad \frac{3\ln^2 2}{176\pi^2}\ \text{(stampato)}.$$
A un loop il minimo è a $x = 1$ esattamente: $b_{\min} = \Lambda^2$; $W_{\min}/\Lambda^4 = -11/(96\pi^2)$.

## 7. La configurazione classica migliorata: il vuoto "spaghetti"

Il tachione condensa in un reticolo bidimensionale di tubi di flusso (Ambjørn–Olesen; rassegna Bordag eq. (146)–(152)). La densità del condensato
$$\rho(z) = e^{-2\pi(\mathrm{Im}\,z)^2/\mathrm{Im}\,\tau}\,|\vartheta_3(\pi z\,|\,\tau)|^2$$
è esattamente periodica sotto le traslazioni magnetiche (verificato). Il parametro di Abrikosov $\beta = \langle\rho^2\rangle/\langle\rho\rangle^2$, calcolato per quadratura diretta dalla definizione: reticolo quadrato $\beta = 1.180340$, triangolare $\beta = 1.159595$ (ottimale; la mutazione obliqua è peggiore). L'energia classica diventa $E_{\rm cl} = \tfrac{H^2}{2}(1 - 1/\beta) < H^2/2$: il reticolo **abbassa** l'energia classica del campo costante.

Riquantizzazione (approssimazione della letteratura: sommare il potenziale a un loop del fondo omogeneo):
$$V(H) = \left(1-\frac1\beta\right)\frac{H^2}{2} + \frac{11}{48\pi^2}(gH)^2\left(\ln\frac{gH}{\mu^2} - \frac12\right),\qquad \ln\frac{gH_{\min}}{\mu^2} = -\left(1-\frac1\beta\right)\frac{24\pi^2}{11g^2},$$
e il minimo è più profondo di Savvidy del fattore $\exp\left(\frac{48\pi^2}{11 g^2\beta}\right)$ — esponenzialmente più profondo a accoppiamento debole.

Nota di convenzione (documentata, non risolta): la eq. (152) della rassegna quota $\vartheta_3(0, q_0) \simeq 1.2713$ per $q_0 = e^{-\pi}$; nella convenzione standard $\vartheta_3(0, e^{-\pi}) = 1.086434$, mentre $1.2713 = \vartheta_3(0, e^{-2})$. La loro eliminazione dei parametri non è ricostruibile dal testo aperto; i nostri numeri sono derivati, non trascritti.

## 8. Verifica del testo primario: Preparata 1986 equazione per equazione

Modulo `preparata_1986.py` (11 test; issue #56, PR #57). Ogni affermazione verificabile del paper è codificata dal suo enunciato testuale più chiaro e ricontrollata con oracoli SymPy, in incrocio con le primitive verificate delle Sezioni 2–6.

1. **Spettro (1.2).** $E^2 = p_z^2 + gH(2n+1) - 2gH\,S_z$: coincide con lo spettro verificato della Sezione 2 (controllo numerico su quattro settori); il settore instabile è $n=0$, $S_z=1$ con $p_z^2 < gH$, come da testo.
2. **Risultato centrale (5.1).** $\Delta E = -\frac{11}{48\pi^2}g^2H^2\ln(\Lambda^2/gH) + O(g^2H^2)$. Il segno **non è trascritto**: è derivato — solo il segno negativo produce il minimo locale con $a, b > 0$ delle (5.2)/(5.3) (il segno positivo dà un massimo, in contraddizione con l'abstract). Il coefficiente del logaritmo è identico al $C = 11/(48\pi^2)$ verificato alla Sezione 2.
3. **Minimo (5.2)/(5.3), derivato.** Dalla (5.1) codificata: $gH^* = a\Lambda^2$ con $a = e^{-1/2}$ e $\Delta E^* = -b\Lambda^4$ con $b = \frac{11}{96\pi^2 e}$ — entrambi finiti per $g \to 0$. Con il termine sotto-leading simbolico $\kappa g^2H^2$: $a = e^{-1/2 - 48\pi^2\kappa/11}$, ancora finito — l'instabilità "essenziale" è robusta, come asserito nel paper.
4. **Algebra AF della Sezione 1.** La (1.1) di Savvidy risolta per $1/g^2$ dà esattamente il running a un loop ($b_0 = 22/3$; residuo nullo). Inserendo il running AF nella variazionale del 1985 (1.3) [$gH^* = \Lambda^2 e^{-12\pi^2/11g^2}$] l'esponente è esattamente $-\ln(\Lambda/\Lambda_{\rm QCD})$, quindi $gH^* = \Lambda\cdot\Lambda_{\rm QCD}$: campo divergente e gap $\sim\Lambda^2$, come affermato.
5. **Distinzione di scala (onestà di convenzione).** La $\Lambda$ di Preparata è il **cut-off UV**; la $\Lambda$ della Sezione 6 è la **scala RG**. I due enunciati ($gH^* = a\Lambda_{\rm UV}^2$ e $b_{\min} = \Lambda_{\rm RG}^2$) non vanno confusi: sono compatibili ma non identici.
6. **Frontiera dichiarata.** La macchina variazionale completa (funzionali d'onda gaussiani, appendici B–F, il modello $\lambda\phi^4$ dell'appendice F) non è riprodotta: è un possibile modulo successivo. Le affermazioni sopra coprono la catena logica portante: spettro → (5.1) → minimo → divergenza del gap → sopravvivenza della conclusione.

**Verdetto:** l'algebra verificabile di Preparata 1986 è internamente consistente e coerente con il coefficiente a un loop verificato indipendentemente; la tesi dell'instabilità essenziale segue dalla (5.1) codificata.

## 9. Risultati dichiarati

1. Il potenziale di Savvidy a un loop è confermato in ogni coefficiente: $C = 11/(48\pi^2)$, minimo a $b = \Lambda^2$, $V_{\min} = -11\Lambda^4/(96\pi^2)$.
2. Il vuoto di Savvidy è **instabile**: tachione per ogni $b>0$, $\mathrm{Im}\,V \neq 0$ a uno e due loop. Il testo primario di Preparata è ora verificato a livello di equazioni (Sezione 8): l'algebra portante è consistente, il minimo $gH^* = e^{-1/2}\Lambda^2$ con $\Delta E^* = -\frac{11}{96\pi^2 e}\Lambda^4$ è derivato dalla sua (5.1), e l'instabilità essenziale ($a$, $b$ finiti per $g \to 0$) ne segue. Frontiera dichiarata: la macchina variazionale completa delle appendici B–F.
3. Il due loop pubblicato (Bordag–Skalozub) contiene tre errori (F1, F2, F3) con i valori corretti esatti sopra; la correzione F3 cambia il coefficiente del termine a due loop e l'esponente del minimo.
4. Il logaritmo a un loop è interamente trasverso-paramagnetico ($2 \times 11/6$); longitudinale e fantasma si cancellano.
5. La configurazione classica più vicina al vuoto quantistico è il reticolo triangolare di tubi di flusso ($\beta = 1.159595$), con energia classica ridotta del fattore $(1 - 1/\beta) = 0.1376$ e minimo esponenzialmente più profondo; la riquantizzazione a un loop attorno ad esso è calcolata sopra.
6. Il potenziale in funzione di $b/\Lambda^2$ è la Sezione 6, con il termine a due loop corretto.

## Riferimenti

[1] G. Preparata, Nuovo Cim. A **96** (1986) 366. DOI: 10.1007/BF02833896 (abstract; testo paywalled).
[2] G.K. Savvidy, Phys. Lett. B **71** (1977) 133; riesposizione aperta: EPJC **80** (2020) 165, arXiv:1910.00654.
[3] N.K. Nielsen, P. Olesen, Nucl. Phys. B **144** (1978) 376.
[4] N.K. Nielsen, P. Olesen, Phys. Lett. B **79** (1978) 304 (linee di vortice elettriche).
[5] J. Ambjørn, P. Olesen, Nucl. Phys. B **170** (1980) 60 e 265 (condensato "spaghetti").
[6] M. Bordag, V. Skalozub, EPJC **82** (2022) 390, arXiv:2112.01043 (aperto).
[7] M. Bordag, Symmetry **15** (2023) 1137 (aperto).
[8] P. Cea, arXiv:2311.14791 (aperto).
[9] Repository: vantasnerdan/substrate-framework — tracker #44, issue #45–#50 e #56, PR #51–#55 e #57; moduli `chromomagnetic_background.py`, `chromomagnetic_two_loop.py`, `chromomagnetic_sectors.py`, `spaghetti_vacuum.py`, `preparata_1986.py` con 63 test; testo primario fissato in `proposals/P229-preparata-qcd-vacuum-audit/sources/`.

*Fine del rapporto.*
