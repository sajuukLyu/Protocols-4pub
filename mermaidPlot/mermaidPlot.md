{
  "theme": "default",
  "themeVariables": {
    "fontSize": "18px",
    "edgeLabelBackground": "#ffffff"
  }
}

```mermaid
flowchart TD;
classDef default stroke-width: 2px ;
classDef qc fill: #ffffff ;
	A("Raw data (.fastq.gz)") -->|fastqc| a{{"Quality Control"}};
	A -->|trimmomatic| B("Trim (.fastq.gz)");
	B -->|fastqc| a;
	B -->|STAR| C("Map (.bam)");
	C -->|RseQC| D{{"Infer experiment"}};
	C -->|featureCount| E("Count (.txt)");
	E --> F("Down stream analysis");
class a,D qc
```