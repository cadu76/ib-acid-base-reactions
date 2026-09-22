// quiz.js - Handles the interactive IB Exam Checkpoint Quiz

const quizQuestions = [
  {
    id: 1,
    question: "A student titrates a 25.0 cm³ sample of 0.10 mol dm⁻³ CH₃COOH (a weak acid, pKₐ = 4.76) with 0.10 mol dm⁻³ NaOH (a strong base). Which of the following indicators is most appropriate to detect the equivalence point of this titration?",
    options: [
      "Methyl Orange (transition range: pH 3.1 - 4.4)",
      "Phenolphthalein (transition range: pH 8.2 - 10.0)",
      "Bromothymol Blue (transition range: pH 6.0 - 7.6)",
      "None of the indicators will show a color change at the equivalence point."
    ],
    answer: 1, // Index 1 is Phenolphthalein
    explanation: "At the equivalence point of a weak acid-strong base titration, a basic salt (sodium acetate, CH₃COONa) is formed. The acetate anion undergoes hydrolysis: CH₃COO⁻ + H₂O ⇌ CH₃COOH + OH⁻, generating hydroxide ions and resulting in an equivalence pH > 7 (specifically around pH 8.73). Phenolphthalein is the most appropriate indicator because its transition range (pH 8.2 - 10.0) encompasses this basic equivalence pH."
  },
  {
    id: 2,
    question: "During the titration of a weak acid with a strong base, the pH at the half-equivalence point is measured to be 4.76. What is the acid dissociation constant (Kₐ) of the weak acid?",
    options: [
      "1.74 × 10⁻⁵ mol dm⁻³",
      "5.75 × 10⁻¹⁰ mol dm⁻³",
      "4.76 mol dm⁻³",
      "1.00 × 10⁻⁷ mol dm⁻³"
    ],
    answer: 0, // Index 0 is 1.74e-5
    explanation: "At the half-equivalence point, exactly 50% of the weak acid (HA) has been neutralized to form its conjugate base (A⁻). Under these buffer conditions, [HA] = [A⁻]. Substituting this equality into the Henderson-Hasselbalch equation (pH = pKₐ + log([A⁻]/[HA])) results in pH = pKₐ. Therefore, pKₐ = 4.76. The dissociation constant Kₐ = 10⁻ᵖᴷᵃ = 10⁻⁴.⁷⁶ = 1.74 × 10⁻⁵ mol dm⁻³."
  },
  {
    id: 3,
    question: "A student measures the electrical conductivity of the solution during the titration of a strong acid (HCl) with a strong base (NaOH). Which statement correctly describes the change in conductivity up to the equivalence point?",
    options: [
      "It increases because more ions are added to the solution.",
      "It remains constant because it is a neutralization reaction.",
      "It decreases because highly mobile H⁺ ions are replaced by less mobile Na⁺ ions.",
      "It drops to zero because water is a non-electrolyte."
    ],
    answer: 2, // Index 2
    explanation: "Protons (H⁺) have exceptionally high molar ionic conductivities due to their proton-hopping mobility mechanism. As NaOH is added, H⁺ reacts with OH⁻ to form neutral water, and is replaced by Na⁺ spectator ions. Because Na⁺ has much lower mobility than H⁺, the overall conductivity of the solution decreases steadily, reaching a minimum at the equivalence point."
  },
  {
    id: 4,
    question: "Why are indicators generally not suitable for determining the equivalence point of a titration between a weak acid (e.g. CH₃COOH) and a weak base (e.g. NH₃)?",
    options: [
      "The equivalence point occurs at a pH of exactly 7, where indicators do not work.",
      "The reaction is not neutral, so no salt is formed.",
      "There is no sharp pH change (inflection jump) near the equivalence point.",
      "Weak acids and weak bases destroy indicator molecules."
    ],
    answer: 2, // Index 2
    explanation: "In a weak acid-weak base titration, the titration curve has a very flat, shallow slope with no sharp vertical pH jump at the equivalence point. Indicators require a sudden pH transition of 2-3 units over a fraction of a milliliter to show a rapid color change. Without this jump, indicator color transitions are gradual and inaccurate."
  },
  {
    id: 5,
    question: "Consider the molecular view of the beaker. If you observe roughly equal amounts of intact acid molecules (HA) and conjugate base anions (A⁻) floating in solution, what is the chemical state of this titration?",
    options: [
      "The initial state before any base was added.",
      "The half-equivalence point of a weak acid titration, forming a buffer.",
      "The equivalence point, where all acid has reacted.",
      "The excess base region after the equivalence point."
    ],
    answer: 1, // Index 1
    explanation: "When [HA] ≈ [A⁻] (intact acid and conjugate base are equal), the system is at the half-equivalence point. This region exhibits maximum buffer capacity, where the solution resists pH changes upon minor additions of acid or base."
  }
];

// Quiz State
let currentQuestionIndex = 0;
let score = 0;
let userAnswers = []; // stores indices of user answers

// UI elements inside the quiz tab
const quizContent = document.getElementById('quizContent');
const quizProgressFill = document.getElementById('quizProgressFill');
const quizProgressText = document.getElementById('quizProgressText');

function initQuiz() {
  currentQuestionIndex = 0;
  score = 0;
  userAnswers = [];
  renderQuestion();
}

function renderQuestion() {
  if (currentQuestionIndex >= quizQuestions.length) {
    renderQuizSummary();
    return;
  }
  
  const q = quizQuestions[currentQuestionIndex];
  
  // Progress Bar
  const progPct = (currentQuestionIndex / quizQuestions.length) * 100;
  quizProgressFill.style.width = `${progPct}%`;
  quizProgressText.textContent = `Question ${currentQuestionIndex + 1} of ${quizQuestions.length}`;
  
  // Build Question Layout
  let optionsHtml = '';
  q.options.forEach((opt, idx) => {
    const letter = String.fromCharCode(65 + idx); // A, B, C, D
    optionsHtml += `
      <button class="opt-btn" onclick="submitAnswer(${idx})">
        <span class="opt-letter">${letter}</span>
        <span>${opt}</span>
      </button>
    `;
  });
  
  quizContent.innerHTML = `
    <div class="question-card">
      <div class="q-header">
        <span class="q-label">Question ${currentQuestionIndex + 1}</span>
        <div class="q-text">${q.question}</div>
      </div>
      
      <div class="options-list">
        ${optionsHtml}
      </div>
      
      <div id="feedbackBox"></div>
    </div>
  `;
}

function submitAnswer(selectedIdx) {
  // Lock further clicks for this question
  document.querySelectorAll('.opt-btn').forEach(btn => {
    btn.classList.add('locked');
    btn.removeAttribute('onclick');
  });
  
  const q = quizQuestions[currentQuestionIndex];
  const buttons = document.querySelectorAll('.opt-btn');
  const feedbackBox = document.getElementById('feedbackBox');
  
  userAnswers.push(selectedIdx);
  const isCorrect = (selectedIdx === q.answer);
  if (isCorrect) score++;
  
  // Color the selected button and correct button
  buttons[selectedIdx].classList.add(isCorrect ? 'correct' : 'incorrect');
  if (!isCorrect) {
    buttons[q.answer].classList.add('correct');
  }
  
  // Render explanation box
  feedbackBox.innerHTML = `
    <div class="quiz-explanation ${isCorrect ? 'success' : 'danger'}">
      <div class="quiz-exp-title ${isCorrect ? 'success' : 'danger'}">
        <i class="fa-solid ${isCorrect ? 'fa-circle-check' : 'fa-circle-xmark'}"></i>
        ${isCorrect ? 'Correct Answer!' : 'Incorrect'}
      </div>
      <div class="quiz-exp-text">${q.explanation}</div>
      
      <div class="quiz-actions">
        <button class="btn-primary" onclick="nextQuestion()">
          ${currentQuestionIndex === quizQuestions.length - 1 ? 'Finish' : 'Next Question'} 
          <i class="fa-solid fa-arrow-right"></i>
        </button>
      </div>
    </div>
  `;
}

function nextQuestion() {
  currentQuestionIndex++;
  renderQuestion();
}

function renderQuizSummary() {
  quizProgressFill.style.width = '100%';
  quizProgressText.textContent = `Completed`;
  
  const percentage = (score / quizQuestions.length) * 100;
  
  let feedbackText = '';
  let icon = 'fa-trophy';
  if (percentage === 100) {
    feedbackText = "Outstanding! Perfect score. You have mastered Grade 12 IB DP Acid-Base Chemistry.";
    icon = 'fa-medal';
  } else if (percentage >= 80) {
    feedbackText = "Excellent job! You have a solid grasp of curves, calculations, and indicators.";
  } else if (percentage >= 50) {
    feedbackText = "Good effort! Review the calculations and molecular animations to strengthen your understanding.";
  } else {
    feedbackText = "Keep practicing! Review the Theory & Calculations section and try the titration simulations again.";
    icon = 'fa-circle-info';
  }
  
  quizContent.innerHTML = `
    <div class="quiz-summary animate-fade-in">
      <div class="summary-icon"><i class="fa-solid ${icon} fa-2xl"></i></div>
      <h2 style="font-size: 1.5rem; margin-top: 0.5rem;">Quiz Results</h2>
      <div class="summary-score">${score} / ${quizQuestions.length}</div>
      <div class="summary-feedback">${feedbackText}</div>
      <div class="summary-desc">You got ${percentage}% of questions correct. Use the Sandbox to reinforce concepts like buffer actions and salt hydrolysis.</div>
      
      <button class="btn-primary" onclick="initQuiz()">
        <i class="fa-solid fa-rotate-left"></i> Retake Quiz
      </button>
    </div>
  `;
}

// Setup on load
window.addEventListener('load', () => {
  initQuiz();
});
