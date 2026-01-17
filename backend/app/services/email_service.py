import os

# Placeholder for Resend
# import resend

class EmailService:
    def __init__(self):
        self.api_key = os.getenv("RESEND_API_KEY")
        # resend.api_key = self.api_key

    async def send_welcome_email(self, to_email: str, name: str):
        """
        [Placeholder] Sends a welcome email using Resend.
        """
        print(f"Sending welcome email to {name} ({to_email})")
        
        # Mock Logic
        # resend.Emails.send({
        #     "from": "onboarding@adc-genai.com",
        #     "to": to_email,
        #     "subject": "Welcome to ADC-GenAI Platform",
        #     "html": "<p>Welcome!</p>"
        # })
        return {"status": "sent"}

    async def send_alert(self, to_email: str, message: str):
        """
        [Placeholder] Sends an alert email.
        """
        print(f"Sending alert to {to_email}: {message}")
        return {"status": "sent"}

email_service = EmailService()
