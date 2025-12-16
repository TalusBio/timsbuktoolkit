use super::*;

#[test]
fn test_app_initialization() {
    let app = ViewerApp::new_test();

    // Check initial state
    assert!(app.data.elution_groups.is_none());
    assert!(matches!(app.data.indexed_data, IndexedDataState::None));
    assert!(app.ui.table_filter.is_empty());
    assert!(app.pending_commands.is_empty());
}

#[test]
fn test_command_handling_select_elution_group() {
    let mut app = ViewerApp::new_test();

    // Create a mock elution group data would be ideal, but for now we just test the command flow
    // Pushing a SelectElutionGroup command
    app.pending_commands
        .push(AppCommand::SelectElutionGroup(42));

    // Process commands
    app.handle_commands();

    // Check that the UI state was updated
    assert_eq!(app.ui.selected_index, Some(42));

    // Check that a RegenerateChromatogram command was chained (as per logic in handle_commands)
    assert_eq!(app.pending_commands.len(), 1);
    assert!(matches!(
        app.pending_commands[0],
        AppCommand::RegenerateChromatogram
    ));
}

#[test]
fn test_command_handling_update_tolerance() {
    let mut app = ViewerApp::new_test();

    // Set a selection so that update tolerance triggers regeneration
    app.ui.selected_index = Some(1);

    app.pending_commands.push(AppCommand::UpdateTolerance);
    app.handle_commands();

    // Should chain regeneration because we have a selection
    assert_eq!(app.pending_commands.len(), 1);
    assert!(matches!(
        app.pending_commands[0],
        AppCommand::RegenerateChromatogram
    ));
}

#[test]
fn test_command_handling_update_tolerance_no_selection() {
    let mut app = ViewerApp::new_test();

    // No selection
    app.ui.selected_index = None;

    app.pending_commands.push(AppCommand::UpdateTolerance);
    app.handle_commands();

    // Should NOT chain regeneration
    assert!(app.pending_commands.is_empty());
}

#[test]
fn test_session_logging() {
    use std::io::Write;
    use std::sync::{
        Arc,
        Mutex,
    };

    // Shared buffer to capture output
    #[derive(Clone)]
    struct SharedWriter(Arc<Mutex<Vec<u8>>>);

    impl Write for SharedWriter {
        fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
            self.0.lock().unwrap().extend_from_slice(buf);
            Ok(buf.len())
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    let shared_buf = Arc::new(Mutex::new(Vec::new()));
    let writer = SharedWriter(shared_buf.clone());

    let mut app = ViewerApp::new_test();
    app.session_log = Some(Box::new(writer));

    // Push a command
    app.pending_commands
        .push(AppCommand::SelectElutionGroup(123));
    app.handle_commands();

    // Check log content
    let output = String::from_utf8(shared_buf.lock().unwrap().clone()).unwrap();

    // Verify JSON structure matches expected serialization of AppCommand::SelectElutionGroup(123)
    // The exact format depends on serde_json, but it should contain the variant name and value
    assert!(output.contains("SelectElutionGroup"));
    assert!(output.contains("123"));
}

#[test]
fn test_auto_zoom_initial_state() {
    let app = ViewerApp::new_test();
    // Should default to 0
    assert_eq!(app.computed.auto_zoom_frame_counter, 0);
}

#[test]
fn test_data_state_initialization() {
    let app = ViewerApp::new_test();
    // Source should be None initially
    assert!(app.data.elution_groups_source.is_none());
}

#[test]
fn test_persistence_save_load() {
    use eframe::{
        App,
        Storage,
    };
    use std::collections::HashMap;
    use std::sync::Mutex; // Import traits!

    // Mock storage
    #[derive(Default)]
    struct MockStorage {
        data: Mutex<HashMap<String, String>>,
    }

    impl eframe::Storage for MockStorage {
        fn get_string(&self, key: &str) -> Option<String> {
            self.data.lock().unwrap().get(key).cloned()
        }

        fn set_string(&mut self, key: &str, value: String) {
            self.data.lock().unwrap().insert(key.to_string(), value);
        }

        fn flush(&mut self) {}
    }

    // 1. Setup initial app and modify state
    let mut app = ViewerApp::new_test();
    app.ui.table_filter = "test_filter".to_string();
    app.file_loader.elution_groups_path = Some(PathBuf::from("/test/path.json"));
    app.ui.selected_index = Some(5);

    // 2. Save
    let mut storage = MockStorage::default();
    app.save(&mut storage);

    // 3. Create context with populated storage
    // We can't easily mock CreationContext cleanly because it has non-public fields or requires specific setup,
    // but ViewerApp::new uses `cc.storage`.
    // Instead of calling ViewerApp::new (which requires CreationContext), let's extract the logic we want to test
    // or simulate what ViewerApp::new does.

    // Simulate ViewerApp::new loading logic:
    let state_string = storage
        .get_string(eframe::APP_KEY)
        .expect("Should have saved state");
    let state: PersistentState = ron::from_str(&state_string).expect("Should deserialize");

    // 4. Verify restored state
    assert_eq!(state.ui_state.table_filter, "test_filter");
    assert_eq!(
        state.file_loader.elution_groups_path,
        Some(PathBuf::from("/test/path.json"))
    );
    assert_eq!(state.ui_state.selected_index, Some(5));
}
