# Algoritmi
Repository for the MSc course: "Algoritmi numerici per la fisica".

# Indicazioni per il codice

Conviene guardare le funzioni scritte nel file "function_alg.cpp" anziché quelle nei vari file di esercizi. Nel file con tutte le funzioni sono presenti molti più commenti e la scrittura e più ordinata.

# Istruzioni per fare il push da terminale

Queste istruzioni valgono dalla seconda volta in poi, devi prima configurare la cartella.

1. Vai nella cartella della Directory (in questo caso fai `cd Algoritmi`)
2. Aggiungi i file al commit:
    ```git add . ```
3. Crea il commit:
    ```git commit -m "Messaggio del commit"```
4. Fai il push:
    ```git push```


# Impostare il Rule a 80 colonne

Il professore richiede, oltre la corretta indentazione, che il codice si possa leggere in una finestra con 80 colonne. Per impostare il ruler ad 80 su VSCode fai:

1. Preferences -> Settings

2. Cerca nella barra di ricerca `rulers`

3. Clicca su `Edit in settings.json` e aggiungi la linea:
    ```"editor.rulers": [80]```