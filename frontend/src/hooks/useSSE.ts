import { useEffect, useRef, useState, useCallback } from 'react';

export interface TaskProgress {
  task_id: string;
  status: 'processing' | 'completed' | 'failed' | 'pending';
  message?: string;
  current?: number;
  total?: number;
  result?: any;
  error?: string;
}

export function useSSE(taskId: string | null, onProgress: (progress: TaskProgress) => void) {
  const [isConnected, setIsConnected] = useState(false);
  const eventSourceRef = useRef<EventSource | null>(null);

  useEffect(() => {
    if (!taskId) {
      return;
    }

    const eventSource = new EventSource(`http://localhost:8000/tasks/${taskId}/stream`);
    eventSourceRef.current = eventSource;

    eventSource.onopen = () => {
      setIsConnected(true);
    };

    eventSource.addEventListener('connected', (event) => {
      const data = JSON.parse(event.data);
      onProgress(data);
    });

    eventSource.addEventListener('progress', (event) => {
      const data = JSON.parse(event.data);
      onProgress(data);
      
      // Close connection if task is complete or failed
      if (data.status === 'completed' || data.status === 'failed') {
        eventSource.close();
        setIsConnected(false);
      }
    });

    eventSource.onerror = (error) => {
      console.error('SSE error:', error);
      setIsConnected(false);
      eventSource.close();
    };

    return () => {
      eventSource.close();
      setIsConnected(false);
    };
  }, [taskId, onProgress]);

  const close = useCallback(() => {
    if (eventSourceRef.current) {
      eventSourceRef.current.close();
      setIsConnected(false);
    }
  }, []);

  return { isConnected, close };
}

