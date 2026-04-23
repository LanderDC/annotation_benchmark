```{mermaid}
flowchart TD
    A([🚀 Start]) --> B[Read command-line arguments\nproteins, categories, output, model, etc.]
    B --> C[Load proteins JSON file]
    C --> D[Load categories JSON file\nand normalize into a clean list]
    D --> E[Load AI language model\nfrom HuggingFace]

    E --> F{Optional:\nQuantize model?}
    F -- 4-bit or 8-bit --> G[Apply memory compression\nto reduce GPU usage]
    F -- No --> H[Load model at full precision]
    G --> I[Model ready]
    H --> I

    I --> J{Does an output file\nalready exist?}
    J -- Yes --> K[Load existing results\nSkip already-processed proteins]
    J -- No --> L[Start fresh]
    K --> M
    L --> M

    M([🔁 Loop over each protein]) --> N{Already\nprocessed?}
    N -- Yes --> M
    N -- No --> O[Extract protein names\nand virus taxonomy]

    O --> P[Build prompt for AI:\nprotein names + taxonomy +\npredefined categories]
    P --> Q[Send prompt to\nAI language model]
    Q --> R[Receive AI response\nin JSON format]

    R --> S{Response\nparseable?}
    S -- Yes --> T[Validate assigned categories\nagainst predefined list]
    S -- No --> U[⚠️ Log warning\nStore empty categories]

    T --> V{Any categories\nneed repair?}
    V -- Yes --> W[Auto-fix minor name issues\ne.g. missing parentheses]
    V -- No --> X
    W --> X[Remove any invalid categories]

    X --> Y[Save result for this protein]
    U --> Y
    Y --> Z{Checkpoint:\nEvery N proteins}
    Z -- Yes --> AA[💾 Save progress to output file]
    Z -- No --> M
    AA --> M

    M -- All proteins done --> AB[💾 Final save to output JSON]
    AB --> AC{Were any original\nlabels provided?}
    AC -- Yes --> AD[Compute accuracy metrics\nPrecision / Recall / F1]
    AC -- No --> AE
    AD --> AE([✅ Done!\nResults saved to output file])

    style A fill:#4CAF50,color:#fff
    style AE fill:#4CAF50,color:#fff
    style M fill:#2196F3,color:#fff
    style Q fill:#9C27B0,color:#fff
    style AA fill:#FF9800,color:#fff
    style AB fill:#FF9800,color:#fff
    style U fill:#f44336,color:#fff
```