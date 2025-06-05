Google-Firebase-Prompt

Create a Firebase web application called "Health Research Summarizer" that allows users to:

1. Upload a PDF file OR paste text from a scientific health research article.
2. Choose from one of three summary types:
   - Plain-language summary (for the general public)
   - Executive summary (for healthcare or policy professionals)
   - Social media summary (Twitter/LinkedIn-style)
3. Submit the content and receive a summary using a large language model.

🔧 Backend:
- Use Firebase Functions to handle the summarization logic.
- Integrate Vertex AI's PaLM 2 or Gemini API (or allow OpenAI via proxy with API key).
- Handle PDF parsing using a Node.js or Python microservice.
- Store summaries in Firebase Firestore for reuse or display history.

📂 Frontend:
- Use Firebase Hosting + a responsive frontend (React, Vue, or plain HTML).
- Include:
   - File upload for PDF
   - Text area for pasted input
   - Dropdown to select summary type
   - "Generate Summary" button
   - Display output summary cleanly

🛡️ Authentication:
- Add optional Firebase Authentication (Google Sign-In)
- Allow anonymous use but limit daily requests without login

📈 Optional Features:
- Track summary usage per user
- Enable side-by-side comparisons (Gemini vs OpenAI vs Claude)
- Provide download/share options for generated summaries

Use Firebase CLI to deploy functions and hosting. All backend code should be modular and secure, with environment variables stored in Firebase config.

