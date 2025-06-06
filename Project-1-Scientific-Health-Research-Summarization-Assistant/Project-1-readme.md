# 🧬 Project 1: Health Research Summarization Assistant

This project explores the use of large language models (LLMs) to automatically summarize complex scientific health research into formats that are easier to understand and share. It is designed for use by healthcare communicators, policy leaders, researchers, and the public.

---

## 🎯 Objective

The goal is to develop an AI-powered summarization tool that converts dense academic papers from journals like *Nature Medicine*, *JAMA*, or *The Lancet* into:

- ✅ Plain-language summaries for general public understanding  
- ✅ Executive bullet-point summaries for healthcare and policy leaders  
- ✅ Social media summaries for science communication on platforms like Twitter and LinkedIn  

---

## 🛠 Key Components

- **Prompt Templates** for different audiences and formats  
- **Model Comparison**: Outputs from OpenAI (GPT-4 / GPT-3.5), Claude, and Gemini  
- **Zero-shot vs Few-shot Prompting** experiments  
- **Chain-of-Thought vs Direct Prompting** experiments  
- **Streamlit Demo App** for uploading and summarizing research PDFs  

---

## 🚀 Getting Started

### Option A: Streamlit App

To run the Streamlit app locally:

1. **Install dependencies:**
   ```bash
   pip install -r requirements.txt

2. **Run the app:**
    streamlit run app.py

3. **Upload a research paper (PDF) or paste plain text**

4. **Choose one of the following summary types:**
    1. Plain Language Summary
    2. Executive Summary
    3. Social Media Summary

5. **Click "Generate Summary" to view the output**

### Option B: LangChain Script:

1. **Install dependencies:**
   ```bash
   pip install -r requirements.txt

2. **Run the app:**
   python run_pipeline.py

3. **Script Will:**
    1. Load a PDF or pre-supplied text
    2. Format it into a prompt
    3. Use GPT-3.5 (or another LLM) to generate three summary formats

---

## 💡 Use Cases

- **Public Health Outreach:**  
  Turn complex studies into accessible language for websites, flyers, or blog posts.

- **Policy Briefing:**  
  Generate bullet summaries for public health and government decision-makers.

- **Science Communication:**  
  Quickly produce accurate summaries for Twitter, LinkedIn, or newsletters.

- **Education & Teaching:**  
  Help students and trainees understand dense medical research papers.

---

## 📊 Experiments

This repository includes controlled prompt experiments to evaluate:

- **Zero-shot vs Few-shot prompting:**  
  Testing the effectiveness of summaries generated without vs. with example prompts.

- **Chain-of-thought vs Direct prompting:**  
  Comparing outputs when models are instructed to reason step-by-step vs. summarize directly.

👉 See detailed results and examples in the `/experiments/` folder.

---

## 🙏 Acknowledgments & Credits

This project uses and builds on several open-source tools and AI services:

- **[OpenAI](https://platform.openai.com/)** — for providing GPT models used in summarization
- **[LangChain](https://www.langchain.com/)** — for LLM chaining and structured prompting
- **[Streamlit](https://streamlit.io/)** — for the live demo frontend interface
- **[PyPDF2](https://pypi.org/project/PyPDF2/)** — for extracting text from uploaded PDFs
- **[PubMed Central](https://www.ncbi.nlm.nih.gov/pmc/)** — for open-access medical research articles
- **[arXiv](https://arxiv.org/)** — for publicly accessible preprints in health and medicine
 
 ---

## 🛡 License

This project is licensed under the **MIT License**.  
See the [LICENSE](./LICENSE) file for full terms of use.

---

## 📝 Citation

If you use this tool or its results in a publication, project, or public communication campaign, please cite:






