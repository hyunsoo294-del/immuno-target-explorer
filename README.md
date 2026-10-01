# Immuno Target Explorer

The default screen is the TAA analysis (gene, cancer, subtype, immune cell, molecule table). It uses one UCSC Xena Toil TPM cohort for patient-tumor bulk RNA and keeps healthy-blood HPA data in a separate reference view. The previous RNA, protein, TCGA, IHC, and character explorer remains under Legacy explorer.

```bat
run_chatbot.bat
```

Gene search for immuno-oncology targets: RNA (HPA/GTEx vs TCGA), protein IHC, TCGA distribution, and CPTAC.

Public sources only. Research use, not for clinical decisions.

## Local run

```bat
run_chatbot.bat
```

Opens http://127.0.0.1:8502

## Share on Streamlit Community Cloud

1. Push this repo to GitHub (public).
2. Open https://share.streamlit.io and sign in with GitHub.
3. **Create app** then this repository, `app.py`, Deploy.

The live URL looks like `https://<name>.streamlit.app`.
