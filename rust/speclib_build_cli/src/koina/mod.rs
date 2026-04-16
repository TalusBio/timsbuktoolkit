pub mod adapter;
pub mod models;

use models::{
    FragmentModel, FragmentPrediction, KoinaRequest, KoinaResponse, PredictionInput, RtModel,
    RtPrediction,
};

// ── Client ───────────────────────────────────────────────────────────────────

pub struct KoinaClient {
    http: reqwest::Client,
    base_url: String,
    fragment_model: FragmentModel,
    rt_model: RtModel,
}

impl KoinaClient {
    pub fn new(base_url: impl Into<String>, fragment_model: FragmentModel, rt_model: RtModel) -> Self {
        Self {
            http: reqwest::Client::new(),
            base_url: base_url.into(),
            fragment_model,
            rt_model,
        }
    }

    /// POST fragment-intensity predictions for `inputs`.
    ///
    /// The batch is sent as a single request; callers are responsible for
    /// splitting into appropriately sized batches before calling this method.
    pub async fn predict_fragments(
        &self,
        inputs: &[PredictionInput],
    ) -> Result<Vec<FragmentPrediction>, String> {
        let url = format!(
            "{}/{}/infer",
            self.base_url,
            self.fragment_model.model_name()
        );
        let request = adapter::build_fragment_request(inputs, "fragment-0");
        let response = self.send_with_retry(&url, &request).await?;
        adapter::parse_fragment_response(response, inputs.len())
    }

    /// POST RT predictions for `inputs`.
    pub async fn predict_rt(
        &self,
        inputs: &[PredictionInput],
    ) -> Result<Vec<RtPrediction>, String> {
        let url = format!("{}/{}/infer", self.base_url, self.rt_model.model_name());
        let request = adapter::build_rt_request(inputs, "rt-0");
        let response = self.send_with_retry(&url, &request).await?;
        adapter::parse_rt_response(response, inputs.len())
    }

    /// Send a request with up to 3 attempts and exponential back-off.
    ///
    /// 4xx responses are not retried (client error, retrying will not help).
    async fn send_with_retry(
        &self,
        url: &str,
        request: &KoinaRequest,
    ) -> Result<KoinaResponse, String> {
        const MAX_ATTEMPTS: u32 = 3;

        let mut last_err = String::new();

        for attempt in 0..MAX_ATTEMPTS {
            if attempt > 0 {
                let delay = std::time::Duration::from_millis(200 * (1u64 << attempt));
                tokio::time::sleep(delay).await;
            }

            let resp = self
                .http
                .post(url)
                .json(request)
                .send()
                .await
                .map_err(|e| e.to_string())?;

            let status = resp.status();

            // Client errors (4xx) — no point retrying.
            if status.is_client_error() {
                let body = resp.text().await.unwrap_or_default();
                return Err(format!(
                    "Koina request failed with {status} (not retrying): {body}"
                ));
            }

            if status.is_success() {
                return resp
                    .json::<KoinaResponse>()
                    .await
                    .map_err(|e| format!("failed to deserialise Koina response: {e}"));
            }

            // 5xx or other non-success — record and retry.
            last_err = format!("Koina request failed with status {status}");
        }

        Err(format!(
            "Koina request failed after {MAX_ATTEMPTS} attempts: {last_err}"
        ))
    }
}
