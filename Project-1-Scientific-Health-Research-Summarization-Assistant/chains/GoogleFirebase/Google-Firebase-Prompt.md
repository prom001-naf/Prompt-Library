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

Fixes/Errors:

This "503 Service Unavailable" error means that the specific AI model you're trying to use (gemini-2.0-flash) is currently overloaded on Google's end. This is usually a temporary issue with the service itself.

While we can't directly fix an overloaded Google server, we can try switching to a different, potentially less busy, compatible model for your summarization tasks. I'll update your Genkit configuration to use gemini-1.5-flash-latest instead of gemini-2.0-flash. This might help circumvent the current overload.

Restore
I've updated the default model in your Genkit configuration. Hopefully, this alternative model will be available and allow your summaries to be generated. If you still encounter issues, it might be necessary to wait a bit and try again, as these kinds of service availability problems are often transient.

Nafis
can you give me correct prompt here

Gemini
The "503 Service Unavailable" error you saw previously is typically an issue on the AI model provider's side (like Google's servers being overloaded) and often not directly related to the prompt content itself, as long as the prompt is valid. We've already switched to the gemini-1.5-flash-latest model, which might help with availability.

However, I can certainly review and refine the prompt in src/ai/flows/summarize-article.ts to ensure it's as clear and effective as possible for the AI. I'll make the instructions for each summary type more detailed and actionable.

Here are the changes:

303582c7
This revised prompt provides more specific guidance for each summary type, which should help the AI generate more targeted and appropriate summaries. If the "503 Service Unavailable" error persists, it's likely due to temporary issues with the AI model provider's service rather than the prompt itself