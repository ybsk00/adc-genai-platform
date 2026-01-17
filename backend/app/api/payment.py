from fastapi import APIRouter, HTTPException, Request
from pydantic import BaseModel

router = APIRouter()

class CheckoutRequest(BaseModel):
    plan_id: str
    user_email: str

@router.post("/create-checkout")
async def create_checkout(data: CheckoutRequest):
    """
    [Placeholder] Creates a Lemon Squeezy checkout session.
    """
    print(f"Creating checkout for {data.user_email} with plan {data.plan_id}")
    
    # Mock Response
    return {
        "checkout_url": "https://store.lemonsqueezy.com/checkout/buy/...",
        "status": "pending"
    }

@router.post("/webhook")
async def payment_webhook(request: Request):
    """
    [Placeholder] Handles Lemon Squeezy webhooks (subscription_created, etc.)
    """
    payload = await request.json()
    print(f"Received webhook: {payload.get('meta', {}).get('event_name')}")
    
    return {"status": "processed"}
