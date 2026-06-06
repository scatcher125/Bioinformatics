# Attempt at learning LangGraph to create agentic workflows
# Reva S
# 06-Jun-2026
# -----------

# 1. Hello world example
from langgraph.graph import StateGraph, MessagesState, START, END

def mock_llm(state: MessagesState):
    return {"messages": [{"role": "ai", "content": "hello world"}]}

graph = StateGraph(MessagesState)
graph.add_node(mock_llm)
graph.add_edge(START, "mock_llm")
graph.add_edge("mock_llm", END)
graph = graph.compile()

response = graph.invoke({"messages": [{"role": "user", "content": "hi!"}]})
print(response)

# modules
# data
# workflow
# save