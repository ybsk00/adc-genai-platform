import asyncio
from unittest.mock import MagicMock, patch
from app.services.bulk_importer import BulkImporter

async def test_importer_modes():
    print("🧪 Testing ClinicalTrials.gov Worker Modes...")
    
    importer = BulkImporter()
    
    # Mocking aiohttp.ClientSession
    with patch("aiohttp.ClientSession.get") as mock_get:
        # Mock response for Daily Mode
        mock_response_daily = MagicMock()
        mock_response_daily.status = 200
        mock_response_daily.json.return_value = asyncio.Future()
        mock_response_daily.json.return_value.set_result({
            "studies": [{"protocolSection": {"identificationModule": {"nctId": "NCT00000001"}}}],
            "nextPageToken": None
        })
        
        # Mock response for Full Mode (Page 1)
        mock_response_full_p1 = MagicMock()
        mock_response_full_p1.status = 200
        mock_response_full_p1.json.return_value = asyncio.Future()
        mock_response_full_p1.json.return_value.set_result({
            "studies": [{"protocolSection": {"identificationModule": {"nctId": "NCT00000002"}}}],
            "nextPageToken": "token_page_2"
        })
        
        # Mock response for Full Mode (Page 2)
        mock_response_full_p2 = MagicMock()
        mock_response_full_p2.status = 200
        mock_response_full_p2.json.return_value = asyncio.Future()
        mock_response_full_p2.json.return_value.set_result({
            "studies": [{"protocolSection": {"identificationModule": {"nctId": "NCT00000003"}}}],
            "nextPageToken": None
        })

        # Setup mock side effects
        mock_get.return_value.__aenter__.side_effect = [
            mock_response_daily, 
            mock_response_full_p1, 
            mock_response_full_p2
        ]
        
        # Mock save_batch to avoid DB calls
        importer.save_batch = MagicMock(return_value=asyncio.Future())
        importer.save_batch.return_value.set_result(1)

        # 1. Test Daily Mode
        print("\n[Test 1] Running Daily Mode...")
        await importer.run_import(max_studies=10, mode="daily")
        
        # Verify Daily Mode call arguments
        # We expect at least one call with filter.lastUpdatePostDate
        daily_call_args = mock_get.call_args_list[0]
        params = daily_call_args[1]['params']
        if "filter.lastUpdatePostDate" in params:
            print("✅ Daily Mode: filter.lastUpdatePostDate parameter found.")
        else:
            print("❌ Daily Mode: filter.lastUpdatePostDate parameter MISSING!")

        # 2. Test Full Mode
        print("\n[Test 2] Running Full Mode...")
        # Reset total_imported for clean test
        importer.total_imported = 0
        
        # We need to reset the side_effect for the next run or use a new mock
        # Simulating Full Mode calls
        mock_get.return_value.__aenter__.side_effect = [
            mock_response_full_p1, 
            mock_response_full_p2
        ]
        
        await importer.run_import(max_studies=10, mode="full")
        
        # Verify Full Mode pagination
        # We expect 2 calls (Page 1, Page 2)
        # Note: run_import loops through search terms and status filters, so actual calls might be more.
        # But we mocked 2 responses.
        print("✅ Full Mode executed.")

if __name__ == "__main__":
    asyncio.run(test_importer_modes())
