# Public deployment with Streamlit Community Cloud

1. Upload the application source to a GitHub repository.
2. Sign in at https://share.streamlit.io using GitHub and choose Create app.
3. Select the repository, its branch, and `app.py` as the entrypoint.
4. In Advanced settings choose Python 3.11 (the locally tested version).
5. Deploy, and set app sharing to public if it is not already public.

The waveguide page is available at `/Waveguide_Modes` on the deployed app URL.
No API keys, secrets or external data services are required.

Required source: app.py, physics/, visualization/, pages/, requirements.txt,
and .streamlit/config.toml. Tests and README.md may also be included.
The local virtual environment, tmp screenshots, credentials, and .git directory
must not be uploaded. Calculation_II.pdf is a development reference, not a
runtime dependency; the deployment bundle omits it.

Official instructions:
https://docs.streamlit.io/deploy/streamlit-community-cloud/deploy-your-app/deploy
