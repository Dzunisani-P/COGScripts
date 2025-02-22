import asyncio
import aiohttp

ids = ["UP000000625", "UP000002311", "UP000000589", "UP000000803"]

async def check_api(session, proteome_id):
    url = f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=proteome:{proteome_id}"
    try:
        async with session.get(url, timeout=10) as response:  # Set timeout for requests
            if response.status == 200:
                return f"{proteome_id}: PASSED"
            else:
                return f"{proteome_id}: FAILED (Status {response.status})"
    except asyncio.TimeoutError:
        return f"{proteome_id}: FAILED (Timeout)"
    except aiohttp.ClientError as e:
        return f"{proteome_id}: FAILED ({e})"

async def main():
    async with aiohttp.ClientSession() as session:
        tasks = [check_api(session, id) for id in ids]
        results = await asyncio.gather(*tasks)

    for result in results:
        print(result)

    if any("FAILED" in res for res in results):
        print("\nThe API's down Bruvv 🚨")
    else:
        print("\nAPI looks good Boiiii ✅")

# Run the event loop
asyncio.run(main())
