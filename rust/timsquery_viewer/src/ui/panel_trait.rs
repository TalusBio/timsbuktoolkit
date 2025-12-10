//! Panel trait and context for consistent UI architecture

use eframe::egui;

use crate::app::{
    AppCommand,
    ComputedState,
    DataState,
    UiState,
};
use crate::file_loader::FileLoader;

/// Write-only command sink
///
/// Prevents panels from seeing or modifying commands from other panels.
/// Panels can only push commands, ensuring proper encapsulation.
pub struct CommandSink<'a> {
    queue: &'a mut Vec<AppCommand>,
}

impl<'a> CommandSink<'a> {
    pub(crate) fn new(queue: &'a mut Vec<AppCommand>) -> Self {
        Self { queue }
    }

    /// Push a command to be processed by the app
    #[inline]
    pub fn push(&mut self, cmd: AppCommand) {
        self.queue.push(cmd);
    }

    /// Push multiple commands at once
    #[inline]
    pub fn extend(&mut self, cmds: impl IntoIterator<Item = AppCommand>) {
        self.queue.extend(cmds);
    }
}

/// Context passed to all panels during rendering
///
/// Provides read access to application state and write-only access to commands.
pub struct PanelContext<'a> {
    pub data: &'a mut DataState,
    pub ui: &'a mut UiState,
    pub computed: &'a ComputedState,
    pub file_loader: &'a mut FileLoader,
    pub commands: CommandSink<'a>,
}

impl<'a> PanelContext<'a> {
    /// Create a new panel context
    pub fn new(
        data: &'a mut DataState,
        ui: &'a mut UiState,
        computed: &'a ComputedState,
        file_loader: &'a mut FileLoader,
        command_queue: &'a mut Vec<AppCommand>,
    ) -> Self {
        Self {
            data,
            ui,
            computed,
            file_loader,
            commands: CommandSink::new(command_queue),
        }
    }
}

/// Trait that all panels must implement
///
/// Provides a consistent interface for rendering UI panels and handling events.
pub trait Panel {
    /// Render the panel
    ///
    /// Panels should read state from `ctx`, render UI using `ui`,
    /// and emit commands by pushing to `ctx.commands`.
    fn render(&mut self, ui: &mut egui::Ui, ctx: &mut PanelContext);

    /// Title displayed in the tab
    fn title(&self) -> &str;
}
