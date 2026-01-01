import { useState, useRef, useEffect, useCallback, useMemo } from 'react';
import { useMolecularStore } from '../store/molecularStore';
import { molecularAPI, ChatResponse, TaskResponse } from '../services/api';
import { useSSE, TaskProgress } from '../hooks/useSSE';
import toast from 'react-hot-toast';

export default function AIChat() {
  const { structure, updateStructure, currentTaskId, setCurrentTaskId, setTaskProgress, taskProgress, clearTask } = useMolecularStore();
  const [messages, setMessages] = useState<Array<{ role: 'user' | 'assistant'; content: string }>>([]);
  const [input, setInput] = useState('');
  const [loading, setLoading] = useState(false);
  const messagesEndRef = useRef<HTMLDivElement>(null);

  // Handle SSE progress updates
  const handleProgress = useCallback((progress: TaskProgress) => {
    setTaskProgress(progress);
    
    if (progress.status === 'completed' && progress.result) {
      setLoading(false);
      
      // Handle chat response
      if (progress.result.response) {
        setMessages(prev => [...prev, { 
          role: 'assistant', content: progress.result.response 
        }]);
        
        // Update structure if provided
        if (progress.result.structure) {
          updateStructure(progress.result.structure);
        }
      }
      
      // Handle structure generation
      if (progress.result.structure && !progress.result.response) {
        updateStructure(progress.result.structure);
        setMessages(prev => [...prev, {
          role: 'assistant',
          content: `I've created the molecule you requested. It's now displayed in the 3D viewer.`
        }]);
      }
      
      clearTask();
    } else if (progress.status === 'failed') {
      setLoading(false);
      toast.error(progress.error || 'Task failed');
      setMessages(prev => [...prev, {
        role: 'assistant',
        content: `Sorry, I encountered an error: ${progress.error || 'Unknown error'}`
      }]);
      clearTask();
    }
  }, [setTaskProgress, updateStructure, clearTask]);

  useSSE(currentTaskId, handleProgress);

  const scrollToBottom = useCallback(() => {
    messagesEndRef.current?.scrollIntoView({ behavior: 'smooth' });
  }, []);

  useEffect(() => {
    scrollToBottom();
  }, [messages, scrollToBottom]);

  const handleSend = useCallback(async () => {
    if (!input.trim() || loading) return;

    const userMessage = input.trim();
    setInput('');
    setMessages(prev => [...prev, { role: 'user', content: userMessage }]);
    setLoading(true);
    clearTask();

    try {
      // Check if this is a structure generation request
      const isGenerationRequest = /create|generate|make|build|show/i.test(userMessage) && 
                                  !structure;
      
      if (isGenerationRequest) {
        // Try to generate structure (may be instant for common molecules)
        try {
          const genResponse = await molecularAPI.generateStructure(userMessage, false);
          if ('structure' in genResponse && genResponse.structure) {
            updateStructure(genResponse.structure);
            setMessages(prev => [...prev, {
              role: 'assistant',
              content: `I've created the molecule you requested. It's now displayed in the 3D viewer.`
            }]);
            setLoading(false);
            return;
          } else if ('error' in genResponse) {
            // Show error and try async
            console.log('Sync generation failed, trying async:', genResponse.error);
          }
        } catch (genError) {
          console.log('Sync generation error:', genError);
        }
        
        // Try async if sync fails or returns no structure
        try {
          const genResponse = await molecularAPI.generateStructure(userMessage, true);
          if ('task_id' in genResponse) {
            const taskResponse = genResponse as TaskResponse;
            setCurrentTaskId(taskResponse.task_id);
            setTaskProgress({
              task_id: taskResponse.task_id,
              status: 'processing',
              message: 'Searching for molecule...',
            });
            return;
          } else if ('structure' in genResponse && genResponse.structure) {
            updateStructure(genResponse.structure);
            setMessages(prev => [...prev, {
              role: 'assistant',
              content: `I've created the molecule you requested. It's now displayed in the 3D viewer.`
            }]);
            setLoading(false);
            return;
          }
        } catch (e) {
          console.log('Async generation failed:', e);
          setMessages(prev => [...prev, {
            role: 'assistant',
            content: 'Sorry, I could not generate that molecule. Please try:\n- A different name or spelling\n- Providing the SMILES notation directly\n- A simpler molecule name'
          }]);
          setLoading(false);
          toast.error('Could not generate molecule');
          return;
        }
      }

      // Regular chat flow - try sync first for fast responses
      try {
        const response = await molecularAPI.chat(userMessage, structure || undefined, false);
        if ('response' in response) {
          const chatResponse = response as ChatResponse;
          setMessages(prev => [...prev, { role: 'assistant', content: chatResponse.response }]);
          if (chatResponse.structure) {
            updateStructure(chatResponse.structure);
          }
          setLoading(false);
          return;
        }
      } catch (e) {
        // Fall through to async if sync fails
      }
      
      // Use async for complex AI queries
      const response = await molecularAPI.chat(userMessage, structure || undefined, true);
      
      if ('task_id' in response) {
        const taskResponse = response as TaskResponse;
        setCurrentTaskId(taskResponse.task_id);
        setTaskProgress({
          task_id: taskResponse.task_id,
          status: 'processing',
          message: 'Processing your message...',
        });
      } else {
        // Synchronous response (fallback)
        const chatResponse = response as ChatResponse;
        setMessages(prev => [...prev, { role: 'assistant', content: chatResponse.response }]);
        if (chatResponse.structure) {
          updateStructure(chatResponse.structure);
        }
        setLoading(false);
      }
    } catch (error) {
      console.error('Chat error:', error);
      setMessages(prev => [...prev, {
        role: 'assistant',
        content: 'Sorry, I encountered an error. Please try again.'
      }]);
      setLoading(false);
      toast.error('Failed to process message');
    }
  }, [input, loading, structure, updateStructure, clearTask, setCurrentTaskId, setTaskProgress]);

  const handleKeyPress = useCallback((e: React.KeyboardEvent) => {
    if (e.key === 'Enter' && !e.shiftKey) {
      e.preventDefault();
      handleSend();
    }
  }, [handleSend]);

  const quickPrompts = useMemo(() => [
    'Create water',
    'Create benzene',
    'Create aspirin',
    'Create caffeine',
  ], []);

  const handleQuickPrompt = useCallback((prompt: string) => {
    setInput(prompt);
  }, []);

  return (
    <div className="flex flex-col h-full bg-transparent">
      {/* Modern Header */}
      <div className="p-4 border-b border-slate-200/50 bg-white/40 backdrop-blur-sm">
        <h3 className="text-sm font-semibold text-slate-900">AI Assistant</h3>
        <p className="text-xs text-slate-500 mt-0.5">Ask me to create or modify molecules</p>
      </div>

      {/* Messages Area */}
      <div className="flex-1 overflow-y-auto p-4 space-y-3 scroll-smooth">
        {messages.length === 0 && (
          <div className="space-y-3">
            <p className="text-xs font-medium text-slate-500 uppercase tracking-wide">Quick Start</p>
            <div className="grid grid-cols-1 gap-2">
              {quickPrompts.map((prompt, idx) => (
                <button
                  key={idx}
                  onClick={() => handleQuickPrompt(prompt)}
                  className="text-left px-3 py-2 bg-white/60 hover:bg-white/80 border border-slate-200/50 rounded-lg text-sm text-slate-700 hover:text-slate-900 transition-all hover:shadow-sm active:scale-[0.98]"
                >
                  {prompt}
                </button>
              ))}
            </div>
          </div>
        )}

        {messages.map((msg, idx) => (
          <div
            key={idx}
            className={`flex ${msg.role === 'user' ? 'justify-end' : 'justify-start'} animate-in fade-in slide-in-from-bottom-2 duration-200`}
          >
            <div
              className={`max-w-[80%] rounded-lg px-3 py-2 text-sm ${
                msg.role === 'user'
                  ? 'bg-gradient-to-r from-blue-500 to-blue-600 text-white shadow-sm'
                  : 'bg-white/80 border border-slate-200/50 text-slate-700 shadow-sm'
              }`}
            >
              <p className="whitespace-pre-wrap leading-relaxed">{msg.content}</p>
            </div>
          </div>
        ))}

        {(loading || taskProgress) && (
          <div className="flex justify-start animate-in fade-in">
            <div className="bg-white/80 border border-slate-200/50 rounded-lg px-3 py-2 max-w-[80%] shadow-sm">
              {taskProgress && taskProgress.message ? (
                <div className="space-y-2">
                  <p className="text-sm text-slate-700">{taskProgress.message}</p>
                  {taskProgress.current !== undefined && taskProgress.total !== undefined && (
                    <>
                      <div className="w-full bg-slate-100 rounded-full h-1.5 overflow-hidden">
                        <div 
                          className="bg-gradient-to-r from-blue-500 to-blue-600 h-full rounded-full transition-all duration-300"
                          style={{ 
                            width: `${(taskProgress.current / taskProgress.total) * 100}%` 
                          }}
                        ></div>
                      </div>
                      <p className="text-xs text-slate-500">
                        {taskProgress.current} / {taskProgress.total}
                      </p>
                    </>
                  )}
                </div>
              ) : (
                <div className="flex space-x-1.5">
                  <div className="w-2 h-2 bg-blue-500 rounded-full animate-bounce" style={{ animationDelay: '0ms' }}></div>
                  <div className="w-2 h-2 bg-blue-500 rounded-full animate-bounce" style={{ animationDelay: '150ms' }}></div>
                  <div className="w-2 h-2 bg-blue-500 rounded-full animate-bounce" style={{ animationDelay: '300ms' }}></div>
                </div>
              )}
            </div>
          </div>
        )}

        <div ref={messagesEndRef} />
      </div>

      {/* Input Area */}
      <div className="p-4 border-t border-slate-200/50 bg-white/40 backdrop-blur-sm">
        <div className="flex gap-2">
          <input
            type="text"
            value={input}
            onChange={(e) => setInput(e.target.value)}
            onKeyPress={handleKeyPress}
            placeholder="Ask me to create or modify a molecule..."
            className="flex-1 bg-white/80 border border-slate-200/50 text-slate-900 rounded-lg px-4 py-2.5 text-sm focus:outline-none focus:ring-2 focus:ring-blue-500 focus:border-transparent transition-all"
            disabled={loading}
          />
          <button
            onClick={handleSend}
            disabled={loading || !input.trim()}
            className="bg-gradient-to-r from-blue-500 to-blue-600 hover:from-blue-600 hover:to-blue-700 disabled:from-slate-300 disabled:to-slate-400 disabled:cursor-not-allowed text-white px-6 py-2.5 rounded-lg text-sm font-medium shadow-sm hover:shadow-md transition-all active:scale-95"
          >
            Send
          </button>
        </div>
      </div>
    </div>
  );
}
