import pandas as pd
import numpy as np
import tensorflow as tf
from sklearn.preprocessing import MultiLabelBinarizer
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dense, Dropout, Conv1D, MaxPooling1D, Flatten, Input, Embedding, LSTM, Concatenate
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.utils import to_categorical
from sklearn.model_selection import train_test_split
from sklearn.metrics import f1_score, precision_score, recall_score, accuracy_score, confusion_matrix, classification_report, matthews_corrcoef
from Bio import SeqIO
import os
import joblib
from sklearn.metrics import f1_score
import seaborn as sns
import matplotlib.pyplot as plt
tf.config.run_functions_eagerly(True)
# Ścieżki do plików Fasta
chromosome_files = {
    "chr1": "/home/paweldyngosz/Dokumenty/Praca_Dyplomowa/Arabidopsis_thaliana.TAIR10.dna.chromosome.1.fa",
    "chr2": "/home/paweldyngosz/Dokumenty/Praca_Dyplomowa/Arabidopsis_thaliana.TAIR10.dna.chromosome.2.fa",
    "chr3": "/home/paweldyngosz/Dokumenty/Praca_Dyplomowa/Arabidopsis_thaliana.TAIR10.dna.chromosome.3.fa",
    "chr4": "/home/paweldyngosz/Dokumenty/Praca_Dyplomowa/Arabidopsis_thaliana.TAIR10.dna.chromosome.4.fa",
    "chr5": "/home/paweldyngosz/Dokumenty/Praca_Dyplomowa/Arabidopsis_thaliana.TAIR10.dna.chromosome.5.fa"
}

# Funkcja wyodrębniania sekwencji z pliku FASTA przy użyciu współrzędnych
def extract_sequence(chromosome_file, start, end):
    print(f"Extracting sequence from {chromosome_file}, Start: {start}, End: {end}")
    with open(chromosome_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            # Ensure start and end are within bounds
            start = max(0, start)
            end = min(len(record.seq), end)
            print(f"Adjusted Start: {start}, Adjusted End: {end}")
            return str(record.seq[start:end])  
    return None  # If no record is found

# Załaduj adnotacje genów
annotation_file = "merged_file.csv"
data = pd.read_csv(annotation_file, sep="\t", encoding="utf-8")
print(f"Data columns: {data.columns}")
data = data[~data["seq_name"].str.lower().str.contains("chrm")]  # Exclude unsupported chromosomes

# Dodaj kolumnę dla sekwencji
promoter_sequences = []
for idx, row in data.iterrows():
    chromosome = row["seq_name"].lower()
    start = max(0, row["start"] - 20)
    end = row["end"]

    # Wyodrębnij sekwencję z odpowiedniego pliku
    chromosome_file = chromosome_files.get(chromosome)
    if chromosome_file:
        try:
            # Użycie funkcji, aby wyodrębnić sekwencję
            sequence = extract_sequence(chromosome_file, start, end)
        except Exception as e:
            print(f"Error extracting sequence for {chromosome} ({start}-{end}): {e}")
            sequence = None
    else:
        print(f"No chromosome file for {chromosome}")
        sequence = None

    promoter_sequences.append(sequence)

# Dodawanie sekwencji do dataframe
data["sequence"] = promoter_sequences

# Sprawdzanie poprawności sekwencji
print(f"Rows with valid sequences: {data['sequence'].notnull().sum()}/{len(data)}")
data = data[data['sequence'].notnull()]  # Keep only rows with valid sequences
print(f"Rows with valid sequences: {len(data)}")
print(data[['seq_name', 'start', 'end', 'sequence']].head())

# Tworzenie końcowej DataFrame
final_df = data[['sequence', 'gene_id', 'GO_ID', 'feature']].copy()
print(final_df.head())
# kodowanie terminów GO
mlb = MultiLabelBinarizer()
to_binary = mlb.fit_transform(final_df['GO_ID'])
print(f"Length of to_binary: {len(to_binary.tolist())}")
print(f"Number of rows in final_df: {len(final_df)}")
print(f"type(to_binary): {type(to_binary)}")
print(f"Shape of to_binary: {to_binary.shape}")
print(to_binary[:5]) 

final_df.copy()

final_df['GO_ID_binary'] = [list(row) for row in to_binary]
final_df = final_df[final_df.loc[:, 'sequence'].notnull()]

# Kodowanie sekwencji za pomocą one-hot encoding
def one_hot_encode_sequence(seq, max_length=20):
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    one_hot = np.zeros((max_length, 4))
    for i, nucleotide in enumerate(seq[:max_length]):
        if nucleotide in mapping:
            one_hot[i, mapping[nucleotide]] = 1
    return one_hot



final_df['sequence_encoded'] = final_df['sequence'].apply(lambda x: one_hot_encode_sequence(x))

# Przygotowanie danych dla modelu
X_seq = np.array(final_df['sequence_encoded'].tolist())
Y = np.array(final_df['GO_ID_binary'].tolist())

print("X_seq shape:", X_seq.shape)
print("Y shape:", Y.shape)

# Train-test split
X_seq_train, X_seq_test, Y_train, Y_test = train_test_split(X_seq, Y, test_size=0.2, random_state=42)

# budowa modelu
input_seq = Input(shape=(X_seq.shape[1], X_seq.shape[2]), name="sequence_input")
x = Conv1D(filters=128, kernel_size=5, activation='relu')(input_seq)
x = MaxPooling1D(pool_size=2)(x)
x = Dropout(0.3)(x)
print(f"Shape before LSTM: {x.shape}")
x = LSTM(128, return_sequences=True)(x)
x = LSTM(64, return_sequences=False)(x)

x = Dense(128, activation='relu')(x)
x = Dropout(0.3)(x)
output = Dense(Y.shape[1], activation='sigmoid', name="output_layer")(x)

model = Model(inputs=[input_seq], outputs=[output])
model.compile(optimizer=Adam(learning_rate=0.001), loss='binary_crossentropy', metrics=['accuracy'])

# Trening modelu
model.fit(X_seq_train, Y_train, validation_split=0.2, epochs=25, batch_size=32, verbose=1)

# Ocena modelu
loss, accuracy = model.evaluate(X_seq_test, Y_test, verbose=0)
print(f"Test Loss: {loss:.4f}, Test Accuracy: {accuracy:.4f}")


# Uzyskje prognozy dla danych testowych
Y_pred = model.predict(X_seq_test)
Y_pred_binary = (Y_pred > 0.5).astype(int)  # Convert probabilities to binary values

# Oblicz confusion matrix
cm = confusion_matrix(np.argmax(Y_test, axis=1), np.argmax(Y_pred_binary, axis=1))

# Wykres
plt.figure(figsize=(10, 7))
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=np.unique(Y_test), yticklabels=np.unique(Y_test))
plt.xlabel('Predicted')
plt.ylabel('True')
plt.title('Confusion Matrix')
plt.show()

# Statystyki
accuracy = accuracy_score(np.argmax(Y_test, axis=1), np.argmax(Y_pred_binary, axis=1))
precision = precision_score(Y_test, Y_pred_binary, average="weighted", zero_division=1)
recall = recall_score(Y_test, Y_pred_binary, average="weighted")
f1 = f1_score(Y_test, Y_pred_binary, average="weighted")
mcc = matthews_corrcoef(np.argmax(Y_test, axis=1), np.argmax(Y_pred_binary, axis=1))
# F1 score
f1 = f1_score(Y_test, Y_pred_binary, average="weighted")  # Use "weighted" for imbalanced classes
print(f"Test F1 Score: {f1:.4f}")


print(f"Test Accuracy: {accuracy:.4f}")
print(f"Test Precision (Weighted): {precision:.4f}")
print(f"Test Recall (Weighted): {recall:.4f}")
print(f"Test F1 Score (Weighted): {f1:.4f}")

print(f"Matthews Correlation Coefficient (MCC): {mcc:.4f}")


model.summary()
# Zapis modelu
model.save("my_model.keras")
with open("go_predictor.pkl", "wb") as f:
    joblib.dump(mlb, f)
